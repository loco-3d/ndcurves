/**
 * \file polynomial.h
 * \brief Definition of a cubic spline.
 * \author Steve T.
 * \version 0.1
 * \date 06/17/2013
 *
 * This file contains definitions for the polynomial struct.
 * It allows the creation and evaluation of natural
 * smooth splines of arbitrary dimension and order
 */

#ifndef _STRUCT_POLYNOMIAL
#define _STRUCT_POLYNOMIAL

#include <algorithm>
#include <functional>
#include <iostream>
#include <stdexcept>

#include "MathDefs.h"
#include "curve_abc.h"

namespace ndcurves {
/// \class polynomial.
/// \brief Represents a polynomial of an arbitrary order defined on the interval
/// \f$[t_{min}, t_{max}]\f$. It follows the equation :<br>
/// \f$ x(t) = a + b(t - t_{min}) + ... + d(t - t_{min})^N \f$<br>
/// where N is the order and \f$ t \in [t_{min}, t_{max}] \f$.
///
template <typename Time = double, typename Numeric = Time, bool Safe = false,
          typename Point = Eigen::Matrix<Numeric, Eigen::Dynamic, 1>,
          typename T_Point =
              std::vector<Point, Eigen::aligned_allocator<Point> > >
struct polynomial : public curve_abc<Time, Numeric, Safe, Point> {
  typedef Point point_t;
  typedef T_Point t_point_t;
  typedef Time time_t;
  typedef Numeric num_t;
  typedef curve_abc<Time, Numeric, Safe, Point> curve_abc_t;
  typedef Eigen::MatrixXd coeff_t;
  typedef Eigen::Ref<coeff_t> coeff_t_ref;
  typedef polynomial<Time, Numeric, Safe, Point, T_Point> polynomial_t;
  typedef typename curve_abc_t::curve_ptr_t curve_ptr_t;

  /* Constructors - destructors */
 public:
  /// \brief Empty constructor. Curve obtained this way can not perform other
  /// class functions.
  ///
  polynomial() : curve_abc_t(), dim_(0), degree_(0), T_min_(0), T_max_(1.0) {}

  /// \brief Constructor.
  /// \param coefficients : a reference to an Eigen matrix where each column is
  /// a coefficient, from the zero order coefficient, up to the highest order.
  /// Spline order is given by the number of the columns -1. \param min  : LOWER
  /// bound on interval definition of the curve. \param max  : UPPER bound on
  /// interval definition of the curve.
  polynomial(const coeff_t& coefficients, const time_t min, const time_t max)
      : curve_abc_t(),
        dim_(coefficients.rows()),
        coefficients_(coefficients),
        degree_(coefficients.cols() - 1),
        T_min_(min),
        T_max_(max) {
    safe_check();
  }

  /// \brief Constructor
  /// \param coefficients : a container containing all coefficients of the
  /// spline, starting
  ///  with the zero order coefficient, up to the highest order. Spline order is
  ///  given by the size of the coefficients.
  /// \param min  : LOWER bound on interval definition of the spline.
  /// \param max  : UPPER bound on interval definition of the spline.
  polynomial(const T_Point& coefficients, const time_t min, const time_t max)
      : curve_abc_t(),
        dim_(coefficients.begin()->size()),
        coefficients_(init_coeffs(coefficients.begin(), coefficients.end())),
        degree_(coefficients_.cols() - 1),
        T_min_(min),
        T_max_(max) {
    safe_check();
  }

  /// \brief Constructor.
  /// \param zeroOrderCoefficient : an iterator pointing to the first element of
  /// a structure containing the coefficients
  ///  it corresponds to the zero degree coefficient.
  /// \param out   : an iterator pointing to the last element of a structure
  /// ofcoefficients. \param min   : LOWER bound on interval definition of the
  /// spline. \param max   : UPPER bound on interval definition of the spline.
  template <typename In>
  polynomial(In zeroOrderCoefficient, In out, const time_t min,
             const time_t max)
      : curve_abc_t(),
        dim_(zeroOrderCoefficient->size()),
        coefficients_(init_coeffs(zeroOrderCoefficient, out)),
        degree_(coefficients_.cols() - 1),
        T_min_(min),
        T_max_(max) {
    safe_check();
  }

  ///
  /// \brief Constructor from boundary condition with C0 : create a polynomial
  /// that connect exactly init and end (order 1) \param init the initial point
  /// of the curve \param end the final point of the curve \param min   : LOWER
  /// bound on interval definition of the spline. \param max   : UPPER bound on
  /// interval definition of the spline.
  ///
  polynomial(const Point& init, const Point& end, const time_t min,
             const time_t max)
      : dim_(init.size()), degree_(1), T_min_(min), T_max_(max) {
    if (T_min_ >= T_max_)
      throw std::invalid_argument("T_min must be strictly lower than T_max");
    if (init.size() != end.size())
      throw std::invalid_argument(
          "init and end points must have the same dimensions.");
    t_point_t coeffs;
    coeffs.push_back(init);
    coeffs.push_back((end - init) / (max - min));
    coefficients_ = init_coeffs(coeffs.begin(), coeffs.end());
    safe_check();
  }

  ///
  /// \brief Constructor from boundary condition with C1 :
  /// create a polynomial that connect exactly init and end and thier first
  /// order derivatives(order 3) \param init the initial point of the curve
  /// \param d_init the initial value of the derivative of the curve
  /// \param end the final point of the curve
  /// \param d_end the final value of the derivative of the curve
  /// \param min   : LOWER bound on interval definition of the spline.
  /// \param max   : UPPER bound on interval definition of the spline.
  ///
  polynomial(const Point& init, const Point& d_init, const Point& end,
             const Point& d_end, const time_t min, const time_t max)
      : dim_(init.size()), degree_(3), T_min_(min), T_max_(max) {
    if (T_min_ >= T_max_)
      throw std::invalid_argument("T_min must be strictly lower than T_max");
    if (init.size() != end.size())
      throw std::invalid_argument(
          "init and end points must have the same dimensions.");
    if (init.size() != d_init.size())
      throw std::invalid_argument(
          "init and d_init points must have the same dimensions.");
    if (init.size() != d_end.size())
      throw std::invalid_argument(
          "init and d_end points must have the same dimensions.");
    /* the coefficients [c0 c1 c2 c3] are found by solving the following system
    of equation (found from the boundary conditions) : [1  0  0   0   ]   [c0]
    [ init ] [1  T  T^2 T^3 ] x [c1] = [ end  ] [0  1  0   0   ]   [c2] [d_init]
    [0  1  2T  3T^2]   [c3]   [d_end ]
    */
    double T = max - min;
    Eigen::Matrix<double, 4, 4> m;
    m << 1., 0, 0, 0, 1., T, T * T, T * T * T, 0, 1., 0, 0, 0, 1., 2. * T,
        3. * T * T;
    Eigen::Matrix<double, 4, 4> m_inv = m.inverse();
    Eigen::Matrix<double, 4, 1> bc;  // boundary condition vector
    coefficients_ = coeff_t::Zero(
        dim_, degree_ + 1);  // init coefficient matrix with the right size
    for (size_t i = 0; i < dim_;
         ++i) {  // for each dimension, solve the boundary condition problem :
      bc[0] = init[i];
      bc[1] = end[i];
      bc[2] = d_init[i];
      bc[3] = d_end[i];
      coefficients_.row(i) = (m_inv * bc).transpose();
    }
    safe_check();
  }

  ///
  /// \brief Constructor from boundary condition with C2 :
  /// create a polynomial that connect exactly init and end and thier first and
  /// second order derivatives(order 5) \param init the initial point of the
  /// curve \param d_init the initial value of the derivative of the curve
  /// \param d_init the initial value of the second derivative of the curve
  /// \param end the final point of the curve
  /// \param d_end the final value of the derivative of the curve
  /// \param d_end the final value of the second derivative of the curve
  /// \param min   : LOWER bound on interval definition of the spline.
  /// \param max   : UPPER bound on interval definition of the spline.
  ///
  polynomial(const Point& init, const Point& d_init, const Point& dd_init,
             const Point& end, const Point& d_end, const Point& dd_end,
             const time_t min, const time_t max)
      : dim_(init.size()), degree_(5), T_min_(min), T_max_(max) {
    if (T_min_ >= T_max_)
      throw std::invalid_argument("T_min must be strictly lower than T_max");
    if (init.size() != end.size())
      throw std::invalid_argument(
          "init and end points must have the same dimensions.");
    if (init.size() != d_init.size())
      throw std::invalid_argument(
          "init and d_init points must have the same dimensions.");
    if (init.size() != d_end.size())
      throw std::invalid_argument(
          "init and d_end points must have the same dimensions.");
    if (init.size() != dd_init.size())
      throw std::invalid_argument(
          "init and dd_init points must have the same dimensions.");
    if (init.size() != dd_end.size())
      throw std::invalid_argument(
          "init and dd_end points must have the same dimensions.");
    /* the coefficients [c0 c1 c2 c3 c4 c5] are found by solving the following
    system of equation (found from the boundary conditions) : [1  0  0   0    0
    0    ]   [c0]   [ init  ] [1  T  T^2 T^3  T^4   T^5  ]   [c1]   [ end   ] [0
    1  0   0    0     0    ]   [c2]   [d_init ] [0  1  2T  3T^2 4T^3  5T^4 ] x
    [c3] = [d_end  ] [0  0  2   0    0     0    ]   [c4]   [dd_init] [0  0  2 6T
    12T^2 20T^3]   [c5]   [dd_end ]
    */
    double T = max - min;
    Eigen::Matrix<double, 6, 6> m;
    m << 1., 0, 0, 0, 0, 0, 1., T, T * T, pow(T, 3), pow(T, 4), pow(T, 5), 0,
        1., 0, 0, 0, 0, 0, 1., 2. * T, 3. * T * T, 4. * pow(T, 3),
        5. * pow(T, 4), 0, 0, 2, 0, 0, 0, 0, 0, 2, 6. * T, 12. * T * T,
        20. * pow(T, 3);
    Eigen::Matrix<double, 6, 6> m_inv = m.inverse();
    Eigen::Matrix<double, 6, 1> bc;  // boundary condition vector
    coefficients_ = coeff_t::Zero(
        dim_, degree_ + 1);  // init coefficient matrix with the right size
    for (size_t i = 0; i < dim_;
         ++i) {  // for each dimension, solve the boundary condition problem :
      bc[0] = init[i];
      bc[1] = end[i];
      bc[2] = d_init[i];
      bc[3] = d_end[i];
      bc[4] = dd_init[i];
      bc[5] = dd_end[i];
      coefficients_.row(i) = (m_inv * bc).transpose();
    }
    safe_check();
  }

  /// \brief Destructor
  virtual ~polynomial() {}

  polynomial(const polynomial& other)
      : dim_(other.dim_),
        coefficients_(other.coefficients_),
        degree_(other.degree_),
        T_min_(other.T_min_),
        T_max_(other.T_max_) {}

  // polynomial& operator=(const polynomial& other);

  /**
   * @brief MinimumJerk Build a polynomial curve connecting p_init to p_final
   * minimizing the time integral of the squared jerk with a zero initial and
   * final velocity and acceleration
   * @param p_init the initial point
   * @param p_final the final point
   * @param t_min initial time
   * @param t_max final time
   * @return the polynomial curve
   */
  static polynomial_t MinimumJerk(const point_t& p_init, const point_t& p_final,
                                  const time_t t_min = 0.,
                                  const time_t t_max = 1.) {
    polynomial_t out =
        polynomial(coeff_t::Zero(p_init.size(), 6), t_min, t_max);
    MinimumJerk(out, p_init, p_final, t_min, t_max);
    return out;
  }

  /**
   * @brief MinimumJerk Build a polynomial curve connecting p_init to p_final
   * minimizing the time integral of the squared jerk with a zero initial and
   * final velocity and acceleration.
   * @param out The output polynomial needs to be of the correct size.
   * @param p_init the initial point
   * @param p_final the final point
   * @param t_min initial time
   * @param t_max final time
   * @return the polynomial curve
   */
  static void MinimumJerk(polynomial_t& out, const point_t& p_init,
                          const point_t& p_final, const time_t t_min = 0.,
                          const time_t t_max = 1.) {
    if (t_min > t_max)
      throw std::invalid_argument(
          "final time should be superior or equal to initial time.");
    const size_t dim(p_init.size());
    if (static_cast<size_t>(p_final.size()) != dim)
      throw std::invalid_argument(
          "Initial and final points must have the same dimension.");
    const double T = t_max - t_min;
    const double T2 = T * T;
    const double T3 = T2 * T;
    const double T4 = T3 * T;
    const double T5 = T4 * T;

    assert(out.coefficients_.cols() == 6);
    assert(out.coefficients_.rows() == static_cast<Eigen::Index>(dim));
    assert(out.dim_ == dim);
    out.coefficients_.fill(0.0);
    out.coefficients_.col(0) = p_init;
    out.coefficients_.col(3) = 10 * (p_final - p_init) / T3;
    out.coefficients_.col(4) = -15 * (p_final - p_init) / T4;
    out.coefficients_.col(5) = 6 * (p_final - p_init) / T5;
    out.degree_ = 5;
    out.T_min_ = t_min;
    out.T_max_ = t_max;
    out.safe_check();
  }

 private:
  void safe_check() {
    if (Safe) {
      if (T_min_ > T_max_) {
        throw std::invalid_argument("Tmin should be inferior to Tmax");
      }
      if (coefficients_.cols() != int(degree_ + 1)) {
        throw std::runtime_error("Spline order and coefficients do not match");
      }
    }
  }

  /* Constructors - destructors */

  /*Operations*/
 public:
  ///  \brief Evaluation of the cubic spline at time t using horner's scheme.
  ///  \param t : time when to evaluate the spline.
  ///  \return \f$x(t)\f$ point corresponding on spline at time t.
  virtual point_t operator()(const time_t t) const {
    check_if_not_empty();
    if ((t < T_min_ || t > T_max_) && Safe) {
      throw std::invalid_argument(
          "error in polynomial : time t to evaluate should be in range [Tmin, "
          "Tmax] of the curve");
    }
    time_t const dt(t - T_min_);
    point_t h = coefficients_.col(degree_);
    for (int i = (int)(degree_ - 1); i >= 0; i--) {
      h = dt * h + coefficients_.col(i);
    }
    return h;
  }

  /**
   * @brief isApprox check if other and *this are approximately equals.
   * Only two curves of the same class can be approximately equals, for
   * comparison between different type of curves see isEquivalent
   * @param other the other curve to check
   * @param prec the precision threshold, default
   * Eigen::NumTraits<Numeric>::dummy_precision()
   * @return true is the two curves are approximately equals
   */
  bool isApprox(
      const polynomial_t& other,
      const Numeric prec = Eigen::NumTraits<Numeric>::dummy_precision()) const {
    return ndcurves::isApprox<num_t>(T_min_, other.min()) &&
           ndcurves::isApprox<num_t>(T_max_, other.max()) &&
           dim_ == other.dim() && degree_ == other.degree() &&
           coefficients_.isApprox(other.coefficients_, prec);
  }

  virtual bool isApprox(
      const curve_abc_t* other,
      const Numeric prec = Eigen::NumTraits<Numeric>::dummy_precision()) const {
    const polynomial_t* other_cast = dynamic_cast<const polynomial_t*>(other);
    if (other_cast)
      return isApprox(*other_cast, prec);
    else
      return false;
  }

  virtual bool operator==(const polynomial_t& other) const {
    return isApprox(other);
  }

  virtual bool operator!=(const polynomial_t& other) const {
    return !(*this == other);
  }

  ///  \brief Evaluation of the derivative of order N of spline at time t.
  ///  \param t : the time when to evaluate the spline.
  ///  \param order : order of derivative.
  ///  \return \f$\frac{d^Nx(t)}{dt^N}\f$ point corresponding on derivative
  ///  spline at time t.
  virtual point_t derivate(const time_t t, const std::size_t order) const {
    check_if_not_empty();
    if ((t < T_min_ || t > T_max_) && Safe) {
      throw std::invalid_argument(
          "error in polynomial : time t to evaluate derivative should be in "
          "range [Tmin, Tmax] of the curve");
    }
    time_t const dt(t - T_min_);
    time_t cdt(1);
    point_t currentPoint_ = point_t::Zero(dim_);
    for (int i = (int)(order); i < (int)(degree_ + 1); ++i, cdt *= dt) {
      currentPoint_ += cdt * coefficients_.col(i) * fact(i, order);
    }
    return currentPoint_;
  }

  polynomial_t compute_derivate(const std::size_t order) const {
    check_if_not_empty();
    if (order == 0) {
      return *this;
    }
    coeff_t coeff_derivated = deriv_coeff(coefficients_);
    polynomial_t deriv(coeff_derivated, T_min_, T_max_);
    return deriv.compute_derivate(order - 1);
  }

  ///  \brief Compute the derived curve at order N.
  ///  \param order : order of derivative.
  ///  \return A pointer to \f$\frac{d^Nx(t)}{dt^N}\f$ derivative order N of the
  ///  curve.
  polynomial_t* compute_derivate_ptr(const std::size_t order) const {
    return new polynomial_t(compute_derivate(order));
  }

  Eigen::MatrixXd coeff() const { return coefficients_; }

  point_t coeffAtDegree(const std::size_t degree) const {
    point_t res;
    if (degree <= degree_) {
      res = coefficients_.col(degree);
    }
    return res;
  }

 private:
  num_t fact(const std::size_t n, const std::size_t order) const {
    num_t res(1);
    for (std::size_t i = 0; i < std::size_t(order); ++i) {
      res *= (num_t)(n - i);
    }
    return res;
  }

  coeff_t deriv_coeff(coeff_t coeff) const {
    if (coeff.cols() == 1)  // only the constant part is left, fill with 0
      return coeff_t::Zero(coeff.rows(), 1);
    coeff_t coeff_derivated(coeff.rows(), coeff.cols() - 1);
    for (std::size_t i = 0; i < std::size_t(coeff_derivated.cols()); i++) {
      coeff_derivated.col(i) = coeff.col(i + 1) * (num_t)(i + 1);
    }
    return coeff_derivated;
  }

  void check_if_not_empty() const {
    if (coefficients_.size() == 0) {
      throw std::runtime_error(
          "Error in polynomial : there is no coefficients set / did you use "
          "empty constructor ?");
    }
  }
  /*Operations*/

 public:
  /*Helpers*/
  /// \brief Get dimension of curve.
  /// \return dimension of curve.
  std::size_t virtual dim() const { return dim_; };
  /// \brief Get the minimum time for which the curve is defined
  /// \return \f$t_{min}\f$ lower bound of time range.
  num_t virtual min() const { return T_min_; }
  /// \brief Get the maximum time for which the curve is defined.
  /// \return \f$t_{max}\f$ upper bound of time range.
  num_t virtual max() const { return T_max_; }
  /// \brief Get the degree of the curve.
  /// \return \f$degree\f$, the degree of the curve.
  virtual std::size_t degree() const { return degree_; }
  /*Helpers*/

  polynomial_t& operator+=(const polynomial_t& p1) {
    assert_operator_compatible(p1);
    if (p1.degree() > degree()) {
      polynomial_t::coeff_t res = p1.coeff();
      res.block(0, 0, coefficients_.rows(), coefficients_.cols()) +=
          coefficients_;
      coefficients_ = res;
      degree_ = p1.degree();
    } else {
      coefficients_.block(0, 0, p1.coeff().rows(), p1.coeff().cols()) +=
          p1.coeff();
    }
    return *this;
  }

  polynomial_t& operator-=(const polynomial_t& p1) {
    assert_operator_compatible(p1);
    if (p1.degree() > degree()) {
      polynomial_t::coeff_t res = -p1.coeff();
      res.block(0, 0, coefficients_.rows(), coefficients_.cols()) +=
          coefficients_;
      coefficients_ = res;
      degree_ = p1.degree();
    } else {
      coefficients_.block(0, 0, p1.coeff().rows(), p1.coeff().cols()) -=
          p1.coeff();
    }
    return *this;
  }

  polynomial_t& operator+=(const polynomial_t::point_t& point) {
    coefficients_.col(0) += point;
    return *this;
  }

  polynomial_t& operator-=(const polynomial_t::point_t& point) {
    coefficients_.col(0) -= point;
    return *this;
  }

  polynomial_t& operator/=(const double d) {
    coefficients_ /= d;
    return *this;
  }

  polynomial_t& operator*=(const double d) {
    coefficients_ *= d;
    return *this;
  }

  ///  \brief Compute the cross product of the current polynomial by another
  ///  polynomial.
  /// The cross product p1Xp2 of 2 polynomials p1 and p2 is defined such that
  /// forall t, p1Xp2(t) = p1(t) X p2(t), with X designing the cross product.
  /// This method of course only makes sense for dimension 3 polynomials.
  ///  \param pOther other polynomial to compute the cross product with.
  ///  \return a new polynomial defining the cross product between this and
  ///  pother
  polynomial_t cross(const polynomial_t& pOther) const {
    assert_operator_compatible(pOther);
    if (dim() != 3)
      throw std::invalid_argument(
          "Can't perform cross product on polynomials with dimensions != 3 ");
    std::size_t new_degree = degree() + pOther.degree();
    coeff_t nCoeffs = Eigen::MatrixXd::Zero(3, new_degree + 1);
    Eigen::Vector3d currentVec;
    Eigen::Vector3d currentVecCrossed;
    for (long i = 0; i < coefficients_.cols(); ++i) {
      currentVec = coefficients_.col(i);
      for (long j = 0; j < pOther.coeff().cols(); ++j) {
        currentVecCrossed = pOther.coeff().col(j);
        nCoeffs.col(i + j) += currentVec.cross(currentVecCrossed);
      }
    }
    // remove last degrees is they are equal to 0
    long final_degree = new_degree;
    while (nCoeffs.col(final_degree).norm() <= ndcurves::MARGIN &&
           final_degree > 0) {
      --final_degree;
    }
    return polynomial_t(nCoeffs.leftCols(final_degree + 1), min(), max());
  }

  ///  \brief Compute the cross product of the current polynomial p by a point
  ///  point.
  /// The cross product pXpoint of is defined such that
  /// forall t, pXpoint(t) = p(t) X point, with X designing the cross product.
  /// This method of course only makes sense for dimension 3 polynomials.
  ///  \param point point to compute the cross product with.
  ///  \return a new polynomial defining the cross product between this and
  ///  point
  polynomial_t cross(const polynomial_t::point_t& point) const {
    if (dim() != 3)
      throw std::invalid_argument(
          "Can't perform cross product on polynomials with dimensions != 3 ");
    coeff_t nCoeffs = coefficients_;
    Eigen::Vector3d currentVec;
    Eigen::Vector3d pointVec = point;
    for (long i = 0; i < coefficients_.cols(); ++i) {
      currentVec = coefficients_.col(i);
      nCoeffs.col(i) = currentVec.cross(pointVec);
    }
    // remove last degrees is they are equal to 0
    long final_degree = degree();
    while (nCoeffs.col(final_degree).norm() <= ndcurves::MARGIN &&
           final_degree > 0) {
      --final_degree;
    }
    return polynomial_t(nCoeffs.leftCols(final_degree + 1), min(), max());
  }

  /*Attributes*/
  std::size_t dim_;       // const
  coeff_t coefficients_;  // const
  std::size_t degree_;    // const
  time_t T_min_, T_max_;  // const
                          /*Attributes*/

 private:
  void assert_operator_compatible(const polynomial_t& other) const {
    if ((fabs(min() - other.min()) > ndcurves::MARGIN) ||
        (fabs(max() - other.max()) > ndcurves::MARGIN) ||
        dim() != other.dim()) {
      throw std::invalid_argument(
          "Can't perform base operation (+ - ) on two polynomials with "
          "different time ranges or different dimensions");
    }
  }

  template <typename In>
  coeff_t init_coeffs(In zeroOrderCoefficient, In highestOrderCoefficient) {
    std::size_t size =
        std::distance(zeroOrderCoefficient, highestOrderCoefficient);
    coeff_t res = coeff_t(dim_, size);
    int i = 0;
    for (In cit = zeroOrderCoefficient; cit != highestOrderCoefficient;
         ++cit, ++i) {
      res.col(i) = *cit;
    }
    return res;
  }

 public:
  // Serialization of the class
  friend class boost::serialization::access;

  template <class Archive>
  void serialize(Archive& ar, const unsigned int version) {
    if (version) {
      // Do something depending on version ?
    }
    ar& BOOST_SERIALIZATION_BASE_OBJECT_NVP(curve_abc_t);
    ar& boost::serialization::make_nvp("dim", dim_);
    ar& boost::serialization::make_nvp("coefficients", coefficients_);
    ar& boost::serialization::make_nvp("dim", dim_);
    ar& boost::serialization::make_nvp("degree", degree_);
    ar& boost::serialization::make_nvp("T_min", T_min_);
    ar& boost::serialization::make_nvp("T_max", T_max_);
  }

};  // class polynomial

template <typename T, typename N, bool S, typename P, typename TP>
polynomial<T, N, S, P, TP> operator+(const polynomial<T, N, S, P, TP>& p1,
                                     const polynomial<T, N, S, P, TP>& p2) {
  polynomial<T, N, S, P, TP> res(p1);
  return res += p2;
}

template <typename T, typename N, bool S, typename P, typename TP>
polynomial<T, N, S, P, TP> operator+(
    const polynomial<T, N, S, P, TP>& p1,
    const typename polynomial<T, N, S, P, TP>::point_t& point) {
  polynomial<T, N, S, P, TP> res(p1);
  return res += point;
}

template <typename T, typename N, bool S, typename P, typename TP>
polynomial<T, N, S, P, TP> operator+(
    const typename polynomial<T, N, S, P, TP>::point_t& point,
    const polynomial<T, N, S, P, TP>& p1) {
  polynomial<T, N, S, P, TP> res(p1);
  return res += point;
}

template <typename T, typename N, bool S, typename P, typename TP>
polynomial<T, N, S, P, TP> operator-(
    const polynomial<T, N, S, P, TP>& p1,
    const typename polynomial<T, N, S, P, TP>::point_t& point) {
  polynomial<T, N, S, P, TP> res(p1);
  return res -= point;
}

template <typename T, typename N, bool S, typename P, typename TP>
polynomial<T, N, S, P, TP> operator-(
    const typename polynomial<T, N, S, P, TP>::point_t& point,
    const polynomial<T, N, S, P, TP>& p1) {
  polynomial<T, N, S, P, TP> res(-p1);
  return res += point;
}

template <typename T, typename N, bool S, typename P, typename TP>
polynomial<T, N, S, P, TP> operator-(const polynomial<T, N, S, P, TP>& p1) {
  typename polynomial<T, N, S, P, TP>::coeff_t res = -p1.coeff();
  return polynomial<T, N, S, P, TP>(res, p1.min(), p1.max());
}

template <typename T, typename N, bool S, typename P, typename TP>
polynomial<T, N, S, P, TP> operator-(const polynomial<T, N, S, P, TP>& p1,
                                     const polynomial<T, N, S, P, TP>& p2) {
  polynomial<T, N, S, P, TP> res(p1);
  return res -= p2;
}

template <typename T, typename N, bool S, typename P, typename TP>
polynomial<T, N, S, P, TP> operator/(const polynomial<T, N, S, P, TP>& p1,
                                     const double k) {
  polynomial<T, N, S, P, TP> res(p1);
  return res /= k;
}

template <typename T, typename N, bool S, typename P, typename TP>
polynomial<T, N, S, P, TP> operator*(const polynomial<T, N, S, P, TP>& p1,
                                     const double k) {
  polynomial<T, N, S, P, TP> res(p1);
  return res *= k;
}

template <typename T, typename N, bool S, typename P, typename TP>
polynomial<T, N, S, P, TP> operator*(const double k,
                                     const polynomial<T, N, S, P, TP>& p1) {
  polynomial<T, N, S, P, TP> res(p1);
  return res *= k;
}

}  // namespace ndcurves

DEFINE_CLASS_TEMPLATE_VERSION(
    SINGLE_ARG(typename Time, typename Numeric, bool Safe, typename Point,
               typename T_Point),
    SINGLE_ARG(ndcurves::polynomial<Time, Numeric, Safe, Point, T_Point>))

#endif  //_STRUCT_POLYNOMIAL
