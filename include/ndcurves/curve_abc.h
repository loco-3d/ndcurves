/**
 * \file curve_abc.h
 * \brief interface for a Curve of arbitrary dimension.
 * \author Steve T.
 * \version 0.1
 * \date 06/17/2013
 *
 * Interface for a curve
 */

#ifndef _STRUCT_CURVE_ABC
#define _STRUCT_CURVE_ABC

#include <functional>
#include <memory>

#include "MathDefs.h"
#include "serialization/archive.hpp"
#include "serialization/eigen-matrix.hpp"
#include "serialization/registeration.hpp"

namespace ndcurves {

template <typename T>
bool isApprox(const T a, const T b, const T eps = 1e-6) {
  return fabs(a - b) < eps;
}

/// \struct curve_abc.
/// \brief Represents a curve of dimension Dim.
/// If value of parameter Safe is false, no verification is made on the
/// evaluation of the curve.
template <typename Time = double, typename Numeric = Time, bool Safe = false,
          typename Point = Eigen::Matrix<Numeric, Eigen::Dynamic, 1>,
          typename Point_derivate = Point>
struct curve_abc : public serialization::Serializable {
  typedef Point point_t;
  typedef Point_derivate point_derivate_t;
  typedef Time time_t;
  typedef Numeric num_t;
  typedef curve_abc<Time, Numeric, Safe, point_t, point_derivate_t>
      curve_t;  // parent class
  typedef curve_abc<Time, Numeric, Safe, point_derivate_t>
      curve_derivate_t;  // parent class
  typedef std::shared_ptr<curve_t> curve_ptr_t;

  /* Constructors - destructors */
 public:
  /// \brief Constructor.
  curve_abc() {}

  /// \brief Destructor.
  virtual ~curve_abc() {}
  /* Constructors - destructors */

  /*Operations*/
  ///  \brief Evaluation of the cubic spline at time t.
  ///  \param t : time when to evaluate the spine
  ///  \return \f$x(t)\f$, point corresponding on curve at time t.
  virtual point_t operator()(const time_t t) const = 0;

  ///  \brief Compute the derived curve at order N.
  ///  \param order : order of derivative.
  ///  \return A pointer to \f$\frac{d^Nx(t)}{dt^N}\f$ derivative order N of the
  ///  curve.
  virtual curve_derivate_t* compute_derivate_ptr(
      const std::size_t order) const = 0;

  /// \brief Evaluate the derivative of order N of curve at time t.
  /// \param t : time when to evaluate the spline.
  /// \param order : order of derivative.
  /// \return \f$\frac{d^Nx(t)}{dt^N}\f$, point corresponding on derivative
  /// curve of order N at time t.
  virtual point_derivate_t derivate(const time_t t,
                                    const std::size_t order) const = 0;

  /**
   * @brief isEquivalent check if other and *this are approximately equal by
   * values, given a precision threshold. This test is done by discretizing both
   * curves and evaluating them and their derivatives.
   * @param other the other curve to check
   * @param order the order up to which the derivatives of the curves are
   * checked for equality
   * @param prec the precision threshold, default
   * Eigen::NumTraits<Numeric>::dummy_precision()
   * @return true if the two curves are approximately equal
   */
  bool isEquivalent(
      const curve_t* other,
      const Numeric prec = Eigen::NumTraits<Numeric>::dummy_precision(),
      const size_t order = 5) const {
    bool equal = ndcurves::isApprox<num_t>(min(), other->min()) &&
                 ndcurves::isApprox<num_t>(max(), other->max()) &&
                 (dim() == other->dim());
    if (!equal) {
      return false;
    }
    time_t inc =
        (max() - min()) / 10.;  // FIXME : define this step somewhere ??
    // check the value along the two curves
    time_t t = min();
    while (t <= max()) {
      if (!(*this)(t).isApprox(other->operator()(t), prec)) {
        return false;
      }
      t += inc;
    }
    //  check if the derivatives are equal
    for (size_t n = 1; n <= order; ++n) {
      t = min();
      while (t <= max()) {
        if (!derivate(t, n).isApprox(other->derivate(t, n), prec)) {
          return false;
        }
        t += inc;
      }
    }
    return true;
  }

  /**
   * @brief isApprox check if other and *this are approximately equal given a
   * precision threshold Only two curves of the same class can be approximately
   * equal, for comparison between different type of curves see isEquivalent.
   * @param other the other curve to check
   * @param prec the precision threshold, default
   * Eigen::NumTraits<Numeric>::dummy_precision()
   * @return true if the two curves are approximately equal
   */
  virtual bool isApprox(
      const curve_t* other,
      const Numeric prec =
          Eigen::NumTraits<Numeric>::dummy_precision()) const = 0;

  /*Operations*/

  /*Helpers*/
  /// \brief Get dimension of curve.
  /// \return dimension of curve.
  virtual std::size_t dim() const = 0;
  /// \brief Get the minimum time for which the curve is defined.
  /// \return \f$t_{min}\f$, lower bound of time range.
  virtual time_t min() const = 0;
  /// \brief Get the maximum time for which the curve is defined.
  /// \return \f$t_{max}\f$, upper bound of time range.
  virtual time_t max() const = 0;
  /// \brief Get the degree of the curve.
  /// \return \f$degree\f$, the degree of the curve.
  virtual std::size_t degree() const = 0;

  std::pair<time_t, time_t> timeRange() { return std::make_pair(min(), max()); }
  /*Helpers*/

  // Serialization of the class
  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version) {
    serialization::register_types<Archive>(ar, version);
    if (version) {
      // Do something depending on version ?
    }
  }
};
BOOST_SERIALIZATION_ASSUME_ABSTRACT(curve_abc)
}  // namespace ndcurves

DEFINE_CLASS_TEMPLATE_VERSION(
    SINGLE_ARG(typename Time, typename Numeric, bool Safe, typename Point,
               typename Point_derivate),
    SINGLE_ARG(ndcurves::curve_abc<Time, Numeric, Safe, Point, Point_derivate>))
#endif  //_STRUCT_CURVE_ABC
