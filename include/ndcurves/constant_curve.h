/**
 * \file constant_curve.h
 * \brief class allowing to create a constant_curve curve.
 * \author Pierre Fernbach
 * \version 0.4
 * \date 29/04/2020
 */

#ifndef _CLASS_CONSTANTCURVE
#define _CLASS_CONSTANTCURVE

#include "curve_abc.h"

namespace ndcurves {
/// \class constant_curve.
/// \brief Represents a constant_curve curve, always returning the same value
/// and a null derivative
///
template <typename Time = double, typename Numeric = Time, bool Safe = false,
          typename Point = Eigen::Matrix<Numeric, Eigen::Dynamic, 1>,
          typename Point_derivate = Point>
struct constant_curve
    : public curve_abc<Time, Numeric, Safe, Point, Point_derivate> {
  typedef Point point_t;
  typedef Point_derivate point_derivate_t;
  typedef Time time_t;
  typedef Numeric num_t;
  typedef constant_curve<Time, Numeric, Safe, Point, Point_derivate>
      constant_curve_t;
  typedef constant_curve<Time, Numeric, Safe, Point_derivate> curve_derivate_t;
  typedef curve_abc<Time, Numeric, Safe, point_t, Point_derivate>
      curve_abc_t;  // parent class

  /* Constructors - destructors */
 public:
  /// \brief Empty constructor. Curve obtained this way can not perform other
  /// class functions.
  ///
  constant_curve() : T_min_(0), T_max_(0), dim_(0) {}

  /// \brief Constructor..
  /// \param value   : The constant value
  /// \param T_min   : lower bound of the time interval
  /// \param T_max   : upper bound of the time interval
  ///
  constant_curve(const Point& value, const time_t T_min = 0.,
                 const time_t T_max = std::numeric_limits<time_t>::max())
      : value_(value), T_min_(T_min), T_max_(T_max), dim_(value.size()) {
    if (Safe && T_min_ > T_max_) {
      throw std::invalid_argument(
          "can't create constant curve: min bound is higher than max bound");
    }
  }

  /// \brief Copy constructor
  /// \param other
  constant_curve(const constant_curve_t& other)
      : value_(other.value_),
        T_min_(other.T_min_),
        T_max_(other.T_max_),
        dim_(other.dim_) {}

  /// \brief Destructor.
  virtual ~constant_curve() {}
  /* Constructors - destructors */

  /*Operations*/
  ///  \brief Evaluation of the cubic spline at time t.
  ///  \param t : time when to evaluate the spine
  ///  \return \f$x(t)\f$, point corresponding on curve at time t.
  virtual point_t operator()(const time_t t) const {
    if (Safe && (t < T_min_ || t > T_max_)) {
      throw std::invalid_argument(
          "error in constant curve : time t to evaluate should be in range "
          "[Tmin, Tmax] of the curve");
    }
    return value_;
  }

  ///  \brief Compute the derived curve at order N.
  ///  Computes the derivative order N, \f$\frac{d^Nx(t)}{dt^N}\f$ of bezier
  ///  curve of parametric equation x(t). \param order : order of derivative.
  ///  \return \f$\frac{d^Nx(t)}{dt^N}\f$ derivative order N of the curve.
  curve_derivate_t compute_derivate() const {
    size_t derivate_size;
    if (point_derivate_t::RowsAtCompileTime == Eigen::Dynamic) {
      derivate_size = dim_;
    } else {
      derivate_size = point_derivate_t::RowsAtCompileTime;
    }
    point_derivate_t value(point_derivate_t::Zero(derivate_size));
    return curve_derivate_t(value, T_min_, T_max_);
  }

  ///  \brief Compute the derived curve at order N.
  ///  \param order : order of derivative.
  ///  \return A pointer to \f$\frac{d^Nx(t)}{dt^N}\f$ derivative order N of the
  ///  curve.
  virtual curve_derivate_t* compute_derivate_ptr(const std::size_t) const {
    return new curve_derivate_t(compute_derivate());
  }

  /// \brief Evaluate the derivative of order N of curve at time t.
  /// \param t : time when to evaluate the spline.
  /// \param order : order of derivative.
  /// \return \f$\frac{d^Nx(t)}{dt^N}\f$, point corresponding on derivative
  /// curve of order N at time t.
  virtual point_derivate_t derivate(const time_t t, const std::size_t) const {
    if (Safe && (t < T_min_ || t > T_max_)) {
      throw std::invalid_argument(
          "error in constant curve : time t to derivate should be in range "
          "[Tmin, Tmax] of the curve");
    }
    size_t derivate_size;
    if (point_derivate_t::RowsAtCompileTime == Eigen::Dynamic) {
      derivate_size = dim_;
    } else {
      derivate_size = point_derivate_t::RowsAtCompileTime;
    }
    return point_derivate_t::Zero(derivate_size);
  }

  /**
   * @brief isApprox check if other and *this are approximately equals given a
   * precision threshold Only two curves of the same class can be approximately
   * equals, for comparison between different type of curves see isEquivalent.
   * @param other the other curve to check
   * @param prec the precision threshold, default
   * Eigen::NumTraits<Numeric>::dummy_precision()
   * @return true is the two curves are approximately equals
   */
  virtual bool isApprox(
      const constant_curve_t& other,
      const Numeric prec = Eigen::NumTraits<Numeric>::dummy_precision()) const {
    return ndcurves::isApprox<num_t>(T_min_, other.min()) &&
           ndcurves::isApprox<num_t>(T_max_, other.max()) &&
           dim_ == other.dim() && value_.isApprox(other.value_, prec);
  }

  virtual bool isApprox(
      const curve_abc_t* other,
      const Numeric prec = Eigen::NumTraits<Numeric>::dummy_precision()) const {
    const constant_curve_t* other_cast =
        dynamic_cast<const constant_curve_t*>(other);
    if (other_cast)
      return isApprox(*other_cast, prec);
    else
      return false;
  }

  virtual bool operator==(const constant_curve_t& other) const {
    return isApprox(other);
  }

  virtual bool operator!=(const constant_curve_t& other) const {
    return !(*this == other);
  }

  /*Helpers*/
  /// \brief Get dimension of curve.
  /// \return dimension of curve.
  std::size_t virtual dim() const { return dim_; }
  /// \brief Get the minimum time for which the curve is defined
  /// \return \f$t_{min}\f$ lower bound of time range.
  num_t virtual min() const { return T_min_; }
  /// \brief Get the maximum time for which the curve is defined.
  /// \return \f$t_{max}\f$ upper bound of time range.
  num_t virtual max() const { return T_max_; }
  /// \brief Get the degree of the curve.
  /// \return \f$degree\f$, the degree of the curve.
  virtual std::size_t degree() const { return 0; }
  /*Helpers*/

  /*Attributes*/
  Point value_;
  time_t T_min_, T_max_;  // const
  std::size_t dim_;       // const
  /*Attributes*/

  // Serialization of the class
  friend class boost::serialization::access;

  template <class Archive>
  void serialize(Archive& ar, const unsigned int version) {
    if (version) {
      // Do something depending on version ?
    }
    ar& BOOST_SERIALIZATION_BASE_OBJECT_NVP(curve_abc_t);
    ar& boost::serialization::make_nvp("value", value_);
    ar& boost::serialization::make_nvp("T_min", T_min_);
    ar& boost::serialization::make_nvp("T_max", T_max_);
    ar& boost::serialization::make_nvp("dim", dim_);
  }

};  // struct constant_curve
}  // namespace ndcurves

DEFINE_CLASS_TEMPLATE_VERSION(
    SINGLE_ARG(typename Time, typename Numeric, bool Safe, typename Point,
               typename Point_derivate),
    SINGLE_ARG(
        ndcurves::constant_curve<Time, Numeric, Safe, Point, Point_derivate>))

#endif  // _CLASS_CONSTANTCURVE
