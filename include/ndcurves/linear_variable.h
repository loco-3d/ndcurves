/**
 * \file linear_variable.h
 * \brief storage for variable points of the form p_i = B_i x + c_i
 * \author Steve T.
 * \version 0.1
 * \date 07/02/2019
 */

#ifndef _CLASS_LINEAR_VARIABLE
#define _CLASS_LINEAR_VARIABLE

#include <math.h>

#include <Eigen/Core>
#include <stdexcept>
#include <vector>

#include "MathDefs.h"
#include "bezier_curve.h"
#include "curve_abc.h"
#include "serialization/archive.hpp"
#include "serialization/eigen-matrix.hpp"

namespace ndcurves {
template <typename Numeric = double, bool Safe = true>
struct linear_variable : public serialization::Serializable {
  typedef Eigen::Matrix<Numeric, Eigen::Dynamic, 1> vector_x_t;
  typedef Eigen::Matrix<Numeric, Eigen::Dynamic, Eigen::Dynamic> matrix_x_t;
  typedef Eigen::Matrix<Numeric, 3, 1> vector_3_t;
  typedef Eigen::Matrix<Numeric, 3, 3> matrix_3_t;
  typedef linear_variable<Numeric> linear_variable_t;

  linear_variable()
      : B_(matrix_x_t::Identity(0, 0)),
        c_(vector_x_t::Zero(0)),
        zero(true) {}  // variable
  linear_variable(const vector_x_t& c)
      : B_(matrix_x_t::Zero(c.size(), c.size())),
        c_(c),
        zero(false) {}  // constant
  linear_variable(const matrix_x_t& B, const vector_x_t& c)
      : B_(B), c_(c), zero(false) {}  // mixed
  linear_variable(const linear_variable_t& other)
      : B_(other.B()),
        c_(other.c()),
        zero(other.isZero()) {}  // copy constructor

  ~linear_variable() {}

  ///  \brief Linear evaluation for vector x.
  ///  \param val : vector to evaluate the linear variable.
  ///  \return Evaluation of linear variable for vector x.
  ///
  vector_x_t operator()(const Eigen::Ref<const vector_x_t>& val) const {
    if (isZero()) return c();
    if (Safe && B().cols() != val.rows())
      throw std::length_error(
          "Cannot evaluate linear variable, variable value does not have the "
          "correct dimension");
    return B() * val + c();
  }

  /// \brief Add another linear variable.
  /// \param w1 : linear variable to add.
  /// \return Linear variable after operation.
  ///
  linear_variable_t& operator+=(const linear_variable_t& w1) {
    if (w1.isZero()) return *this;
    if (isZero()) {
      this->B_ = w1.B_;
      this->c_ = w1.c_;
      zero = w1.isZero();
    } else {
      if (Safe && B().rows() != w1.B().rows())
        throw std::length_error(
            "Cannot add linear variables, variables do not have the same "
            "dimension");
      else if (B().cols() >
               w1.B().cols()) {  // new variables added left for primitive
        B_.block(0, B().cols() - w1.B().cols(), B().rows(), w1.B().cols()) +=
            w1.B();
        c_.tail(w1.c().rows()) += w1.c();
      } else if (B().cols() <
                 w1.B().cols()) {  // new variables added left for primitive
        linear_variable_t opp = w1 + (*this);
        this->B_ = opp.B_;
        this->c_ = opp.c_;
      } else {
        this->B_ += w1.B_;
        this->c_ += w1.c_;
      }
    }
    return *this;
  }

  /// \brief Substract another linear variable.
  /// \param w1 : linear variable to substract.
  /// \return Linear variable after operation.
  ///
  linear_variable_t& operator-=(const linear_variable_t& w1) {
    if (w1.isZero()) return *this;
    if (isZero()) {
      this->B_ = -w1.B_;
      this->c_ = -w1.c_;
      zero = w1.isZero();
    } else {
      if (Safe && B().rows() != w1.B().rows())
        throw std::length_error(
            "Cannot add linear variables, variables do not have the same "
            "dimension");
      else if (B().cols() >
               w1.B().cols()) {  // new variables added left for primitive
        B_.block(0, B().cols() - w1.B().cols(), B().rows(), w1.B().cols()) -=
            w1.B();
        c_.tail(w1.c().rows()) -= w1.c();
      } else if (B().cols() <
                 w1.B().cols()) {  // new variables added left for primitive
        linear_variable_t opp = -w1 + (*this);
        this->B_ = opp.B_;
        this->c_ = opp.c_;
      } else {
        this->B_ -= w1.B_;
        this->c_ -= w1.c_;
      }
    }
    return *this;
  }

  /// \brief Divide by a constant : p_i / d = B_i*x/d + c_i/d.
  /// \param d : constant.
  /// \return Linear variable after operation.
  ///
  linear_variable_t& operator/=(const double d) {
    B_ /= d;
    c_ /= d;
    return *this;
  }

  /// \brief Multiply by a constant : p_i / d = B_i*x*d + c_i*d.
  /// \param d : constant.
  /// \return Linear variable after operation.
  ///
  linear_variable_t& operator*=(const double d) {
    B_ *= d;
    c_ *= d;
    return *this;
  }

  ///  \brief Compute the cross product of the current linear_variable and the
  ///  other.
  /// This method of course only makes sense for dimension 3 curves and
  /// dimension 3 unknown, since otherwise the result is non-linear. It assumes
  /// that a method point_t cross(const point_t&, const point_t&) has been
  /// defined
  ///  \param pOther other polynomial to compute the cross product with.
  ///  \return a new polynomial defining the cross product between this and
  ///  other
  linear_variable_t cross(const linear_variable_t& other) const {
    if (B().rows() != 3)
      throw std::invalid_argument(
          "Can't perform cross product on linear variables with dimensions != "
          "3 ");
    if (B().cols() != 3)
      throw std::invalid_argument(
          "Can't perform cross product on linear variables more than one "
          "unknown ");
    if (isZero() || other.isZero()) return linear_variable_t::Zero(3);
    if ((B().squaredNorm() - B().diagonal().squaredNorm() > MARGIN) ||
        (other.B().squaredNorm() - other.B().diagonal().squaredNorm() > MARGIN))
      throw std::invalid_argument(
          "Can't perform cross product on linear variables if B is not "
          "diagonal ");
    // (B1 x + c1) X (B2 x + c2) = (-c2X B1) x + (bX B2) x + b1Xb2
    typename linear_variable_t::matrix_3_t newB =
        skew<typename linear_variable_t::matrix_3_t,
             typename linear_variable_t::vector_3_t>(-other.c()) *
            B() +
        skew<typename linear_variable_t::matrix_3_t,
             typename linear_variable_t::vector_3_t>(c()) *
            other.B();
    typename linear_variable_t::vector_3_t newC =
        ndcurves::cross(c(), other.c());
    return linear_variable_t(newB, newC);
  }

  /// \brief Get a linear variable equal to zero.
  /// \param dim : Dimension of linear variable.
  /// \return Linear variable equal to zero.
  ///
  static linear_variable_t Zero(size_t dim = 0) {
    return linear_variable_t(matrix_x_t::Zero(dim, dim), vector_x_t::Zero(dim));
  }

  /// \brief Get a linear variable equal to the variable
  /// \param dim : Dimension of linear variable.
  /// \return Linear variable equal to the variable.
  ///
  static linear_variable_t X(size_t dim = 0) {
    return linear_variable_t(matrix_x_t::Identity(dim, dim),
                             vector_x_t::Zero(dim));
  }

  /// \brief Get dimension of linear variable.
  /// \return Dimension of linear variable.
  ///
  std::size_t size() const { return zero ? 0 : std::max(B_.rows(), c_.size()); }

  /// \brief Get norm of linear variable (Norm of B plus norm of C).
  /// \return Norm of linear variable.
  Numeric norm() const { return isZero() ? 0 : (B_.norm() + c_.norm()); }

  /// \brief Check if actual linear variable and other are approximately equal
  /// given a precision threshold. Only two curves of the same class can be
  /// approximately equal, \param prec : the precision threshold, default
  /// Eigen::NumTraits<Numeric>::dummy_precision() \return true if the two
  /// linear variables are approximately equal.
  bool isApprox(
      const linear_variable_t& other,
      const double prec = Eigen::NumTraits<Numeric>::dummy_precision()) const {
    return (*this - other).norm() < prec;
  }

  /// \brief Check if actual linear variable and other are equal.
  /// \param other : the other linear_variable to check.
  /// \return true if the two linear_variable are equals.
  virtual bool operator==(const linear_variable& other) const {
    return this->B_ == other.B_ && this->c_ == other.c_;
  }

  /// \brief Check if actual linear variable and other are different.
  /// \param other : the other linear_variable to check.
  /// \return true if the two linear_variable are different.
  virtual bool operator!=(const linear_variable& other) const {
    return !(*this == other);
  }

  const matrix_x_t& B() const { return B_; }
  const vector_x_t& c() const { return c_; }
  bool isZero() const { return zero; }

  // Serialization of the class
  friend class boost::serialization::access;

  template <class Archive>
  void serialize(Archive& ar, const unsigned int version) {
    if (version) {
      // Do something depending on version ?
    }
    ar& boost::serialization::make_nvp("B_", B_);
    ar& boost::serialization::make_nvp("c_", c_);
    ar& boost::serialization::make_nvp("zero", zero);
  }

  linear_variable& operator=(const linear_variable& other) {
    if (this == &other) {
      return *this;
    }
    // Perform a deep copy here to copy all necessary data.
    // Make sure to handle memory allocation properly.
    // You may need to copy the data contained within the linear_variable.
    this->B_ = other.B_;
    this->c_ = other.c_;
    this->zero = other.zero;
    return *this;
  }

 private:
  matrix_x_t B_;
  vector_x_t c_;
  bool zero;
};

template <typename N, bool S>
inline linear_variable<N, S> operator+(const linear_variable<N, S>& w1,
                                       const linear_variable<N, S>& w2) {
  linear_variable<N, S> res(w1.B(), w1.c());
  return res += w2;
}

template <typename N, bool S>
linear_variable<N, S> operator-(const linear_variable<N, S>& w1,
                                const linear_variable<N, S>& w2) {
  linear_variable<N, S> res(w1.B(), w1.c());
  return res -= w2;
}

template <typename N, bool S>
linear_variable<N, S> operator-(const linear_variable<N, S>& w1) {
  return linear_variable<N, S>(-w1.B(), -w1.c());
}

template <typename N, bool S>
linear_variable<N, S> operator*(const double k,
                                const linear_variable<N, S>& w) {
  linear_variable<N, S> res(w.B(), w.c());
  return res *= k;
}

template <typename N, bool S>
linear_variable<N, S> operator*(const linear_variable<N, S>& w,
                                const double k) {
  linear_variable<N, S> res(w.B(), w.c());
  return res *= k;
}

template <typename N, bool S>
linear_variable<N, S> operator/(const linear_variable<N, S>& w,
                                const double k) {
  linear_variable<N, S> res(w.B(), w.c());
  return res /= k;
}

template <typename BezierFixed, typename BezierLinear, typename X>
BezierFixed evaluateLinear(const BezierLinear& bIn, const X x) {
  typename BezierFixed::t_point_t fixed_wps;
  for (typename BezierLinear::cit_point_t cit = bIn.waypoints().begin();
       cit != bIn.waypoints().end(); ++cit)
    fixed_wps.push_back(cit->operator()(x));
  return BezierFixed(fixed_wps.begin(), fixed_wps.end(), bIn.T_min_,
                     bIn.T_max_);
}

template <typename N, bool S>
std::ostream& operator<<(std::ostream& os, const linear_variable<N, S>& l) {
  return os << "linear_variable: \n \t B:\n"
            << l.B() << "\t c: \n"
            << l.c().transpose();
}

}  // namespace ndcurves

DEFINE_CLASS_TEMPLATE_VERSION(
    SINGLE_ARG(typename Numeric, bool Safe),
    SINGLE_ARG(ndcurves::linear_variable<Numeric, Safe>))
#endif  //_CLASS_LINEAR_VARIABLE
