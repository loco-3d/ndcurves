/**
 * \file quadratic_variable.h
 * \brief storage for variable points of the form p_i = x' A_i x + B_i x + c_i
 * \author Steve T.
 * \version 0.1
 * \date 07/02/2019
 */

#ifndef _CLASS_QUADRATIC_VARIABLE
#define _CLASS_QUADRATIC_VARIABLE

#include <math.h>

#include <Eigen/Core>
#include <stdexcept>
#include <vector>

#include "MathDefs.h"
#include "ndcurves/curve_abc.h"
#include "ndcurves/linear_variable.h"

namespace ndcurves {

template <typename Numeric = double>
struct quadratic_variable {
  typedef Eigen::Matrix<Numeric, Eigen::Dynamic, Eigen::Dynamic> matrix_x_t;
  typedef Eigen::Matrix<Numeric, Eigen::Dynamic, 1> point_t;
  typedef quadratic_variable<Numeric> quadratic_variable_t;

  quadratic_variable() {
    c_ = 0.;
    b_ = point_t::Zero(1);
    A_ = matrix_x_t::Zero(1, 1);
    zero = true;
  }

  quadratic_variable(const matrix_x_t& A, const point_t& b, const Numeric c = 0)
      : c_(c), b_(b), A_(A), zero(false) {
    if (A.cols() != b.rows() || A.cols() != A.rows()) {
      throw std::invalid_argument("The dimensions of A and b are incorrect.");
    }
  }

  quadratic_variable(const point_t& b, const Numeric c = 0)
      : c_(c),
        b_(b),
        A_(matrix_x_t::Zero((int)(b.rows()), (int)(b.rows()))),
        zero(false) {}

  static quadratic_variable_t Zero(size_t /*dim*/ = 0) {
    return quadratic_variable_t();
  }

  // linear evaluation
  Numeric operator()(const Eigen::Ref<const point_t>& val) const {
    if (isZero()) {
      throw std::runtime_error("Not initialized! (isZero)");
    }
    return val.transpose() * A() * val + b().transpose() * val + c();
  }

  quadratic_variable& operator+=(const quadratic_variable& w1) {
    if (w1.isZero()) return *this;
    if (isZero()) {
      this->A_ = w1.A_;
      this->b_ = w1.b_;
      this->c_ = w1.c_;
      zero = false;
    } else {
      this->A_ += w1.A_;
      this->b_ += w1.b_;
      this->c_ += w1.c_;
    }
    return *this;
  }
  quadratic_variable& operator-=(const quadratic_variable& w1) {
    if (w1.isZero()) return *this;
    if (isZero()) {
      this->A_ = -w1.A_;
      this->b_ = -w1.b_;
      this->c_ = -w1.c_;
      zero = false;
    } else {
      this->A_ -= w1.A_;
      this->b_ -= w1.b_;
      this->c_ -= w1.c_;
    }
    return *this;
  }

  quadratic_variable& operator/=(const double d) {
    // handling zero case
    if (!isZero()) {
      this->A_ /= d;
      this->b_ /= d;
      this->c_ /= d;
    }
    return *this;
  }
  quadratic_variable& operator*=(const double d) {
    // handling zero case
    if (!isZero()) {
      this->A_ *= d;
      this->b_ *= d;
      this->c_ *= d;
    }
    return *this;
  }

  const matrix_x_t& A() const {
    if (isZero()) {
      throw std::runtime_error("Not initialized! (isZero)");
    }
    return A_;
  }
  const point_t& b() const {
    if (isZero()) {
      throw std::runtime_error("Not initialized! (isZero)");
    }
    return b_;
  }
  const Numeric c() const {
    if (isZero()) {
      throw std::runtime_error("Not initialized! (isZero)");
    }
    return c_;
  }
  bool isZero() const { return zero; }
  std::size_t size() const {
    return zero ? 0 : std::max(A_.cols(), (std::max(b_.cols(), c_.size())));
  }

 private:
  Numeric c_;
  point_t b_;
  matrix_x_t A_;
  bool zero;
};

/// \brief Transforms a vector into a diagonal matrix
template <typename N>
Eigen::Matrix<N, Eigen::Dynamic, Eigen::Dynamic> to_diagonal(
    const Eigen::Ref<const Eigen::Matrix<N, Eigen::Dynamic, 1> > vec) {
  typedef typename Eigen::Matrix<N, Eigen::Dynamic, Eigen::Dynamic> matrix_t;
  return vec.asDiagonal();
  matrix_t res(matrix_t::Zero(vec.rows(), vec.rows()));
  for (int i = 0; i < vec.rows(); ++i) res(i, i) = vec(i);
  return res;
}

// only works with diagonal linear variables
template <typename N>
inline quadratic_variable<N> operator*(const linear_variable<N>& w1,
                                       const linear_variable<N>& w2) {
  typedef quadratic_variable<N> quad_var_t;
  typedef linear_variable<N> lin_var_t;
  typedef typename quad_var_t::matrix_x_t matrix_x_t;
  typedef typename quad_var_t::point_t point_t;
  typedef typename lin_var_t::vector_x_t point_dim_t;
  point_dim_t ones = point_dim_t::Ones(w1.c().size());
  point_t b1 = w1.B().transpose() * ones, b2 = w2.B().transpose() * ones;
  matrix_x_t B1 = to_diagonal<N>(b1);
  matrix_x_t B2 = to_diagonal<N>(b2);  // b1.array().square()
  // omitting all transposes since A matrices are diagonals
  matrix_x_t A = B1.transpose() * B2;
  point_t b = w1.c().transpose() * w2.B() + w2.c().transpose() * w1.B();
  N c = w1.c().transpose() * w2.c();
  return quad_var_t(A, b, c);
}

template <typename N>
inline quadratic_variable<N> operator+(const quadratic_variable<N>& w1,
                                       const quadratic_variable<N>& w2) {
  quadratic_variable<N> res(w1.A(), w1.b(), w1.c());
  return res += w2;
}

template <typename N>
quadratic_variable<N> operator-(const quadratic_variable<N>& w1,
                                const quadratic_variable<N>& w2) {
  quadratic_variable<N> res(w1.A(), w1.b(), w1.c());
  return res -= w2;
}

template <typename N>
quadratic_variable<N> operator*(const double k,
                                const quadratic_variable<N>& w) {
  quadratic_variable<N> res(w.A(), w.b(), w.c());
  return res *= k;
}

template <typename N>
quadratic_variable<N> operator*(const quadratic_variable<N>& w,
                                const double k) {
  quadratic_variable<N> res(w.A(), w.b(), w.c());
  return res *= k;
}

template <typename N>
quadratic_variable<N> operator/(const quadratic_variable<N>& w,
                                const double k) {
  quadratic_variable<N> res(w.A(), w.b(), w.c());
  return res /= k;
}

}  // namespace ndcurves
#endif  //_CLASS_QUADRATIC_VARIABLE
