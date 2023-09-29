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

#ifndef NDCURVES_POLYNOMIAL_1D_H
#define NDCURVES_POLYNOMIAL_1D_H

#include "ndcurves/polynomial.h"

namespace ndcurves {
/// \class polynomial.
/// \brief Represents a polynomial of an arbitrary order defined on the interval
/// \f$[t_{min}, t_{max}]\f$. It follows the equation :<br>
/// \f$ x(t) = a + b(t - t_{min}) + ... + d(t - t_{min})^N \f$<br>
/// where N is the order and \f$ t \in [t_{min}, t_{max}] \f$.
///
template <typename Time, typename Numeric, bool Safe>
struct polynomial<Time, Numeric, Safe, double, std::vector<double> >
    : public curve_abc<Time, Numeric, Safe, double> {
  typedef double point_t;
  typedef std::vector<double> t_point_t;
  typedef Time time_t;
  typedef Numeric num_t;
  typedef curve_abc<Time, Numeric, Safe, point_t> curve_abc_t;
  typedef Eigen::VectorXd coeff_t;
  typedef Eigen::Ref<coeff_t> coeff_t_ref;
  typedef polynomial<Time, Numeric, Safe, point_t, t_point_t> polynomial_t;
  typedef typename curve_abc_t::curve_ptr_t curve_ptr_t;

  /* Constructors - destructors */
 public:
  /// \brief Empty constructor. Curve obtained this way can not perform other
  /// class functions.
  ///
  polynomial()
      : curve_abc_t(),
        dim_(1),
        t_min_(0),
        t_max_(1),
        degree_(0),
        coefficients_() {}

  /// \brief Constructor.
  /// \param coefficients : a reference to an Eigen matrix where each column is
  /// a coefficient, from the zero order coefficient, up to the highest order.
  /// Spline order is given by the number of the columns -1. \param min  : LOWER
  /// bound on interval definition of the curve. \param max  : UPPER bound on
  /// interval definition of the curve.
  polynomial(const coeff_t& coefficients, const time_t min, const time_t max)
      : curve_abc_t(),
        dim_(1),
        coefficients_(coefficients),
        degree_(coefficients.size() - 1),
        t_min_(min),
        t_max_(max) {
    safe_check();
  }

  /// \brief Constructor
  /// \param coefficients : a container containing all coefficients of the
  /// spline, starting
  ///  with the zero order coefficient, up to the highest order. Spline order is
  ///  given by the size of the coefficients.
  /// \param min  : LOWER bound on interval definition of the spline.
  /// \param max  : UPPER bound on interval definition of the spline.
  polynomial(const t_point_t& coefficients, const time_t min, const time_t max)
      : curve_abc_t(),
        dim_(1),
        coefficients_(Eigen::Map<Eigen::VectorXd>(coefficients.data(),
                                                  coefficients.size())),
        degree_(coefficients_.size() - 1),
        t_min_(min),
        t_max_(max) {
    safe_check();
  }

  /// \brief Constructor.
  /// \param zeroOrderCoefficient : an iterator pointing to the first element of
  /// a structure containing the coefficients
  ///  it corresponds to the zero degree coefficient.
  /// \param out   : an iterator pointing to the last element of a structure
  /// of coefficients. \param min   : LOWER bound on interval definition of the
  /// spline. \param max   : UPPER bound on interval definition of the spline.
  template <typename In>
  polynomial(In zeroOrderCoefficient, In out, const time_t min,
             const time_t max)
      : curve_abc_t(),
        dim_(zeroOrderCoefficient->size()),
        coefficients_(Eigen::Map<Eigen::VectorXd>(
            zeroOrderCoefficient, std::distance(zeroOrderCoefficient, out))),
        degree_(coefficients_.size() - 1),
        t_min_(min),
        t_max_(max) {
    safe_check();
  }

  ///
  /// \brief Constructor from boundary condition with C0 : create a polynomial
  /// that connect exactly init and end (order 1) \param init the initial point
  /// of the curve \param end the final point of the curve \param min   : LOWER
  /// bound on interval definition of the spline. \param max   : UPPER bound on
  /// interval definition of the spline.
  ///
  polynomial(const point_t& init, const point_t& end, const time_t min,
             const time_t max)
      : dim_(1), degree_(1), t_min_(min), t_max_(max) {
    if (t_min_ >= t_max_)
      throw std::invalid_argument("T_min must be strictly lower than T_max");
    coefficients_.setZero(2);
    coefficients_ << init, (end - init) / (max - min);
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
  polynomial(const point_t& init, const point_t& d_init, const point_t& end,
             const point_t& d_end, const time_t min, const time_t max)
      : dim_(1), degree_(3), t_min_(min), t_max_(max) {
    if (t_min_ >= t_max_)
      throw std::invalid_argument("T_min must be strictly lower than T_max");
    /** the coefficients [c0 c1 c2 c3] are found by solving the following system
     * of equation (found from the boundary conditions) :
     * [1  0  0   0   ] [c0]   [ init ]
     * [1  T  T^2 T^3 ] [c1] = [ end  ]
     * [0  1  0   0   ] [c2]   [d_init]
     * [0  1  2T  3T^2] [c3]   [d_end ]
     */
    double T = max - min;
    Eigen::Matrix<double, 4, 4> m;
    m << 1., 0, 0, 0, 1., T, T * T, T * T * T, 0, 1., 0, 0, 0, 1., 2. * T,
        3. * T * T;
    Eigen::Matrix<double, 4, 4> m_inv = m.inverse();
    Eigen::Matrix<double, 4, 1> bc;              // boundary condition vector
    coefficients_ = coeff_t::Zero(degree_ + 1);  // init coefficient size
    for (size_t i = 0; i < dim_;
         ++i) {  // for each dimension, solve the boundary condition problem :
      bc[0] = init;
      bc[1] = end;
      bc[2] = d_init;
      bc[3] = d_end;
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
  polynomial(const point_t& init, const point_t& d_init, const point_t& dd_init,
             const point_t& end, const point_t& d_end, const point_t& dd_end,
             const time_t min, const time_t max)
      : dim_(1), degree_(5), t_min_(min), t_max_(max) {
    if (t_min_ >= t_max_)
      throw std::invalid_argument("T_min must be strictly lower than T_max");
    /**
     * The coefficients [c0 c1 c2 c3 c4 c5] are found by solving the following
     * system of equation (found from the boundary conditions) :
     * [1  0  0   0    0     0    ]   [c0]   [ init  ]
     * [1  T  T^2 T^3  T^4   T^5  ]   [c1]   [ end   ]
     * [0  1  0   0    0     0    ]   [c2]   [d_init ]
     * [0  1  2T  3T^2 4T^3  5T^4 ] x [c3] = [d_end  ]
     * [0  0  2   0    0     0    ]   [c4]   [dd_init]
     * [0  0  2   6T   12T^2 20T^3]   [c5]   [dd_end ]
     */
    double T = max - min;
    Eigen::Matrix<double, 6, 6> m;
    m << 1., 0, 0, 0, 0, 0, 1., T, T * T, pow(T, 3), pow(T, 4), pow(T, 5), 0,
        1., 0, 0, 0, 0, 0, 1., 2. * T, 3. * T * T, 4. * pow(T, 3),
        5. * pow(T, 4), 0, 0, 2, 0, 0, 0, 0, 0, 2, 6. * T, 12. * T * T,
        20. * pow(T, 3);
    Eigen::Matrix<double, 6, 6> m_inv = m.inverse();
    Eigen::Matrix<double, 6, 1> bc;  // boundary condition vector
    coefficients_.setZero(degree_ + 1);
    bc(0) = init;
    bc(1) = end;
    bc(2) = d_init;
    bc(3) = d_end;
    bc(4) = dd_init;
    bc(5) = dd_end;
    coefficients_ = (m_inv * bc).transpose();
    safe_check();
  }

  /// \brief Destructor
  virtual ~polynomial() {}

  polynomial(const polynomial& other)
      : dim_(other.dim_),
        coefficients_(other.coefficients_),
        degree_(other.degree_),
        t_min_(other.t_min_),
        t_max_(other.t_max_) {}

  // polynomial& operator=(const polynomial& other);

  /**
   * @brief This implements a minimum jerk trajectory between two points with
   *        zero velocity and acceleration at the extremities.
   *        This is a real-time safe method, hence no dynamic allocation is
   *        made. This imply that the user already set the size of the
   *        coefficient matrix of the polynomial to the correct size, i.e. 6.
   *
   * @param p_init Initial point from which the trajectory starts.
   * @param p_final Reaching point of the trajectory.
   * @param t_min Starting time point.
   * @param t_max Stop time point.
   */
  void generate_minimum_jerk(const point_t& p_init, const point_t& p_final,
                             const time_t t_min = 0., const time_t t_max = 1.) {
    if (t_min > t_max) {
      throw std::invalid_argument(
          "final time should be superior or equal to initial time.");
    }
    assert(coefficients_.size() == 6);

    const double T = t_max - t_min;
    const double T2 = T * T;
    const double T3 = T2 * T;
    const double T4 = T3 * T;
    const double T5 = T4 * T;

    coefficients_.fill(0.0);
    coefficients_(0) = p_init;
    coefficients_(3) = 10 * (p_final - p_init) / T3;
    coefficients_(4) = -15 * (p_final - p_init) / T4;
    coefficients_(5) = 6 * (p_final - p_init) / T5;

    dim_ = 1;
    degree_ = 5;
    t_min_ = t_min;
    t_max_ = t_max;
    safe_check();
    return;
  }

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
    if (t_min > t_max) {
      throw std::invalid_argument(
          "final time should be superior or equal to initial time.");
    }
    const double T = t_max - t_min;
    const double T2 = T * T;
    const double T3 = T2 * T;
    const double T4 = T3 * T;
    const double T5 = T4 * T;

    coeff_t coeffs = coeff_t::Zero(6);
    coeffs(0) = p_init;
    coeffs(3) = 10 * (p_final - p_init) / T3;
    coeffs(4) = -15 * (p_final - p_init) / T4;
    coeffs(5) = 6 * (p_final - p_init) / T5;
    return polynomial_t(coeffs, t_min, t_max);
  }

 private:
  void safe_check() {
    if (Safe) {
      if (t_min_ > t_max_) {
        throw std::invalid_argument("Tmin should be inferior to Tmax");
      }
      if (coefficients_.size() != int(degree_ + 1)) {
        throw std::runtime_error("Spline order and coefficients do not match");
      }
    }
  }

  /* Constructors - destructors */

  /*Operations*/
 public:
  ///  \brief Evaluation of the cubic spline at time t using horner's scheme.
  ///         Implementation for point_t = double.
  ///  \param t : time when to evaluate the spline.
  ///  \return \f$x(t)\f$ point corresponding on spline at time t.
  point_t operator()(const time_t t) const {
    check_if_not_empty();
    if ((t < t_min_ || t > t_max_) && Safe) {
      throw std::invalid_argument(
          "error in polynomial : time t to evaluate should be in range [Tmin, "
          "Tmax] of the curve");
    }
    point_t out(0.0);
    time_t const dt(t - t_min_);
    out = coefficients_(degree_);
    for (int i = static_cast<int>(degree_ - 1); i >= 0; i--) {
      out = dt * out + coefficients_(i);
    }
    return out;
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
    return ndcurves::isApprox<num_t>(t_min_, other.min()) &&
           ndcurves::isApprox<num_t>(t_max_, other.max()) &&
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
    if ((t < t_min_ || t > t_max_) && Safe) {
      throw std::invalid_argument(
          "error in polynomial : time t to evaluate derivative should be in "
          "range [Tmin, Tmax] of the curve");
    }
    point_t out;
    time_t const dt(t - t_min_);
    time_t cdt(1);
    out = 0.0;
    for (int i = (int)(order); i < (int)(degree_ + 1); ++i, cdt *= dt) {
      out += cdt * coefficients_(i) * fact(i, order);
    }
    return out;
  }

  polynomial_t compute_derivate(const std::size_t order) const {
    check_if_not_empty();
    if (order == 0) {
      return *this;
    }
    coeff_t coeff_derivated = deriv_coeff(coefficients_);
    polynomial_t deriv(coeff_derivated, t_min_, t_max_);
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
      res = coefficients_(degree);
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
  num_t virtual min() const { return t_min_; }
  /// \brief Get the maximum time for which the curve is defined.
  /// \return \f$t_{max}\f$ upper bound of time range.
  num_t virtual max() const { return t_max_; }
  /// \brief Get the degree of the curve.
  /// \return \f$degree\f$, the degree of the curve.
  virtual std::size_t degree() const { return degree_; }
  /*Helpers*/

  polynomial_t& operator+=(const polynomial_t& p1) {
    assert_operator_compatible(p1);
    if (p1.degree() > degree()) {
      polynomial_t::coeff_t res = p1.coeff();
      res.block(0, 0, coefficients_.rows(), coefficients_.size()) +=
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
      res.block(0, 0, coefficients_.rows(), coefficients_.size()) +=
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
    coefficients_(0) += point;
    return *this;
  }

  polynomial_t& operator-=(const polynomial_t::point_t& point) {
    coefficients_(0) -= point;
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
    for (long i = 0; i < coefficients_.size(); ++i) {
      currentVec = coefficients_(i);
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
    for (long i = 0; i < coefficients_.size(); ++i) {
      currentVec = coefficients_(i);
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
  time_t t_min_, t_max_;  // const
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
    ar& boost::serialization::make_nvp("T_min", t_min_);
    ar& boost::serialization::make_nvp("T_max", t_max_);
  }

};  // class polynomial

}  // namespace ndcurves

#endif  // NDCURVES_POLYNOMIAL_1D_H