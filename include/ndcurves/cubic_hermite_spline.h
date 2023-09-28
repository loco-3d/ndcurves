/**
 * \file cubic_hermite_spline.h
 * \brief class allowing to create a cubic hermite spline of any dimension.
 * \author Justin Carpentier <jcarpent@laas.fr> modified by Jason Chemin
 * <jchemin@laas.fr> \date 05/2019
 */

#ifndef _CLASS_CUBICHERMITESPLINE
#define _CLASS_CUBICHERMITESPLINE

#include <boost/serialization/utility.hpp>  // To serialize std::pair
#include <iostream>
#include <stdexcept>
#include <vector>

#include "MathDefs.h"
#include "bezier_curve.h"
#include "curve_abc.h"
#include "piecewise_curve.h"

namespace ndcurves {
/// \class CubicHermiteSpline.
/// \brief Represents a set of cubic hermite splines defining a continuous
/// function \f$p(t)\f$. A hermite cubic spline is a minimal degree polynom
/// interpolating a function in two points \f$P_i\f$ and \f$P_{i+1}\f$ with its
/// tangent \f$m_i\f$ and \f$m_{i+1}\f$.<br> A hermite cubic spline :
/// - crosses each of the waypoint given in its initialization (\f$P_0\f$,
/// \f$P_1\f$,...,\f$P_N\f$).
/// - has its derivatives on \f$P_i\f$ and \f$P_{i+1}\f$ are \f$p'(t_{P_i}) =
/// m_i\f$ and \f$p'(t_{P_{i+1}}) = m_{i+1}\f$.
///
template <typename Time = double, typename Numeric = Time, bool Safe = false,
          typename Point = Eigen::Matrix<Numeric, Eigen::Dynamic, 1> >
struct cubic_hermite_spline : public curve_abc<Time, Numeric, Safe, Point> {
  typedef Point point_t;
  typedef std::pair<Point, Point> pair_point_tangent_t;
  typedef std::vector<pair_point_tangent_t,
                      Eigen::aligned_allocator<pair_point_tangent_t> >
      t_pair_point_tangent_t;
  typedef std::vector<Time> vector_time_t;
  typedef Time time_t;
  typedef Numeric num_t;
  typedef curve_abc<Time, Numeric, Safe, point_t> curve_abc_t;  // parent class
  typedef cubic_hermite_spline<Time, Numeric, Safe, point_t>
      cubic_hermite_spline_t;
  typedef bezier_curve<Time, Numeric, Safe, point_t> bezier_t;
  typedef typename bezier_t::t_point_t t_point_t;
  typedef piecewise_curve<Time, Numeric, Safe, point_t, point_t, bezier_t>
      piecewise_bezier_t;

 public:
  /// \brief Empty constructor. Curve obtained this way can not perform other
  /// class functions.
  ///
  cubic_hermite_spline() : dim_(0), T_min_(0), T_max_(0) {}

  /// \brief Constructor.
  /// \param PairsBegin : an iterator pointing to the first element of a
  /// pair(position, derivative) container. \param PairsEnd   : an iterator
  /// pointing to the last  element of a pair(position, derivative) container.
  /// \param time_control_points : vector containing time for each waypoint.
  ///
  template <typename In>
  cubic_hermite_spline(In PairsBegin, In PairsEnd,
                       const vector_time_t& time_control_points)
      : size_(std::distance(PairsBegin, PairsEnd)), degree_(3) {
    // Check size of pairs container.
    if (Safe && size_ < 1) {
      throw std::length_error(
          "can not create cubic_hermite_spline, number of pairs is inferior to "
          "2.");
    }
    // Set dimension according to size of points
    dim_ = PairsBegin->first.size();
    // Push all pairs in controlPoints
    In it(PairsBegin);
    for (; it != PairsEnd; ++it) {
      if (Safe && (static_cast<size_t>(it->first.size()) != dim_ ||
                   static_cast<size_t>(it->second.size()) != dim_))
        throw std::invalid_argument(
            "All the control points and their derivatives must have the same "
            "dimension.");
      control_points_.push_back(*it);
    }
    // Set time
    setTime(time_control_points);
  }

  cubic_hermite_spline(const cubic_hermite_spline& other)
      : dim_(other.dim_),
        control_points_(other.control_points_),
        time_control_points_(other.time_control_points_),
        duration_splines_(other.duration_splines_),
        T_min_(other.T_min_),
        T_max_(other.T_max_),
        size_(other.size_),
        degree_(other.degree_) {}

  /// \brief Destructor.
  virtual ~cubic_hermite_spline() {}

  /*Operations*/
 public:
  ///  \brief Evaluation of the cubic hermite spline at time t.
  ///  \param t : time when to evaluate the spline.
  ///  \return \f$p(t)\f$ point corresponding on spline at time t.
  ///
  virtual Point operator()(const time_t t) const {
    check_conditions();
    if (Safe & !(T_min_ <= t && t <= T_max_)) {
      throw std::invalid_argument(
          "can't evaluate cubic hermite spline, out of range");
    }
    if (size_ == 1) {
      return control_points_.front().first;
    } else {
      const bezier_t bezier = buildCurrentBezier(t);
      return bezier(t);
    }
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
      const cubic_hermite_spline_t& other,
      const Numeric prec = Eigen::NumTraits<Numeric>::dummy_precision()) const {
    bool equal = ndcurves::isApprox<num_t>(T_min_, other.min()) &&
                 ndcurves::isApprox<num_t>(T_max_, other.max()) &&
                 dim_ == other.dim() && degree_ == other.degree() &&
                 size_ == other.size() &&
                 time_control_points_ == other.time_control_points_ &&
                 duration_splines_ == other.duration_splines_;
    if (!equal) return false;
    for (std::size_t i = 0; i < size_; ++i) {
      if ((!control_points_[i].first.isApprox(other.control_points_[i].first,
                                              prec)) ||
          (!control_points_[i].second.isApprox(other.control_points_[i].second,
                                               prec)))
        return false;
    }
    return true;
  }

  virtual bool isApprox(
      const curve_abc_t* other,
      const Numeric prec = Eigen::NumTraits<Numeric>::dummy_precision()) const {
    const cubic_hermite_spline_t* other_cast =
        dynamic_cast<const cubic_hermite_spline_t*>(other);
    if (other_cast)
      return isApprox(*other_cast, prec);
    else
      return false;
  }

  virtual bool operator==(const cubic_hermite_spline_t& other) const {
    return isApprox(other);
  }

  virtual bool operator!=(const cubic_hermite_spline_t& other) const {
    return !(*this == other);
  }

  ///  \brief Evaluate the derivative of order N of spline at time t.
  ///  \param t : time when to evaluate the spline.
  ///  \param order : order of derivative.
  ///  \return \f$\frac{d^Np(t)}{dt^N}\f$ point corresponding on derivative
  ///  spline of order N at time t.
  ///
  virtual Point derivate(const time_t t, const std::size_t order) const {
    check_conditions();
    if (Safe & !(T_min_ <= t && t <= T_max_)) {
      throw std::invalid_argument(
          "can't derivate cubic hermite spline, out of range");
    }
    if (size_ == 1) {
      return control_points_.front().second;
    } else {
      const bezier_t bezier = buildCurrentBezier(t);
      return bezier.derivate(t, order);
    }
  }

  piecewise_bezier_t compute_derivate(const std::size_t order) const {
    piecewise_bezier_t res;
    for (size_t i = 0; i < size_ - 1; ++i) {
      const bezier_t curve = buildCurrentBezier(time_control_points_[i]);
      res.add_curve(curve.compute_derivate(order));
    }
    return res;
  }

  ///  \brief Compute the derived curve at order N.
  ///  \param order : order of derivative.
  ///  \return A pointer to \f$\frac{d^Nx(t)}{dt^N}\f$ derivative order N of the
  ///  curve.
  piecewise_bezier_t* compute_derivate_ptr(const std::size_t order) const {
    return new piecewise_bezier_t(compute_derivate(order));
  }

  /// \brief Set time of each control point of cubic hermite spline.
  /// Set duration of each spline, Exemple : \f$( 0., 0.5, 0.9, ..., 4.5 )\f$
  /// with values corresponding to times for \f$P_0, P_1, P_2, ..., P_N\f$
  /// respectively.<br> \param time_control_points : Vector containing time for
  /// each control point.
  ///
  void setTime(const vector_time_t& time_control_points) {
    time_control_points_ = time_control_points;
    T_min_ = time_control_points_.front();
    T_max_ = time_control_points_.back();
    if (time_control_points.size() != size()) {
      throw std::length_error(
          "size of time control points should be equal to number of control "
          "points");
    }
    computeDurationSplines();
    if (!checkDurationSplines()) {
      throw std::invalid_argument(
          "time_splines not monotonous, all spline duration should be superior "
          "to 0");
    }
  }

  /// \brief Get vector of pair (positition, derivative) corresponding to
  /// control points. \return vector containing control points.
  ///
  t_pair_point_tangent_t getControlPoints() { return control_points_; }

  /// \brief Get vector of Time corresponding to Time for each control point.
  /// \return vector containing time of each control point.
  ///
  vector_time_t getTime() { return time_control_points_; }

  /// \brief Get number of control points contained in the trajectory.
  /// \return number of control points.
  ///
  std::size_t size() const { return size_; }

  /// \brief Get number of intervals (subsplines) contained in the trajectory.
  /// \return number of intervals (subsplines).
  ///
  std::size_t numIntervals() const { return size() - 1; }

 private:
  /// \brief Get index of the interval (subspline) corresponding to time t for
  /// the interpolation. \param t : time where to look for interval. \return
  /// Index of interval for time t.
  ///
  std::size_t findInterval(const time_t t) const {
    // time before first control point time.
    if (t <= time_control_points_[0]) {
      return 0;
    }
    // time is after last control point time
    if (t >= time_control_points_[size_ - 1]) {
      return size_ - 2;
    }
    std::size_t left_id = 0;
    std::size_t right_id = size_ - 1;
    while (left_id <= right_id) {
      const std::size_t middle_id = left_id + (right_id - left_id) / 2;
      if (time_control_points_.at(middle_id) < t) {
        left_id = middle_id + 1;
      } else if (time_control_points_.at(middle_id) > t) {
        right_id = middle_id - 1;
      } else {
        return middle_id;
      }
    }
    return left_id - 1;
  }

  /**
   * @brief buildCurrentBezier set up the current_bezier_ attribut to represent
   * the curve of the interval that contain t. This bezier is defined by the
   * following control points: p0, p0 + m0/3, p1 - m1/3, p1
   * @param t the time for which the bezier is build
   * @return the bezier curve
   */
  bezier_t buildCurrentBezier(const time_t t) const {
    size_t id_interval = findInterval(t);
    const pair_point_tangent_t pair0 = control_points_.at(id_interval);
    const pair_point_tangent_t pair1 = control_points_.at(id_interval + 1);
    const Time& t0 = time_control_points_[id_interval];
    const Time& t1 = time_control_points_[id_interval + 1];
    t_point_t control_points;
    control_points.reserve(4);
    control_points.push_back(pair0.first);
    control_points.push_back(pair0.first + pair0.second / 3. * (t1 - t0));
    control_points.push_back(pair1.first - pair1.second / 3. * (t1 - t0));
    control_points.push_back(pair1.first);
    return bezier_t(control_points.begin(), control_points.end(), t0, t1);
  }

  /// \brief Check if control points list is not empty and dimension of point
  /// superior to zero.
  ///
  void check_conditions() const {
    if (control_points_.size() == 0) {
      throw std::runtime_error(
          "Error in cubic hermite : there is no control points set / did you "
          "use empty constructor ?");
    } else if (dim_ == 0) {
      throw std::runtime_error(
          "Error in cubic hermite : Dimension of points is zero / did you use "
          "empty constructor ?");
    }
  }

  /// \brief compute duration of each spline.
  /// For N control points with time \f$T_{P_0}, T_{P_1}, T_{P_2}, ...,
  /// T_{P_N}\f$ respectively, Duration of each subspline is : (
  /// T_{P_1}-T_{P_0}, T_{P_2}-T_{P_1}, ..., T_{P_N}-T_{P_{N-1} ).
  ///
  void computeDurationSplines() {
    duration_splines_.clear();
    Time actual_time;
    Time prev_time = *(time_control_points_.begin());
    std::size_t i = 0;
    for (i = 0; i < size() - 1; i++) {
      actual_time = time_control_points_.at(i + 1);
      duration_splines_.push_back(actual_time - prev_time);
      prev_time = actual_time;
    }
  }

  /// \brief Check if duration of each subspline is strictly positive.
  /// \return true if all duration of strictly positive, false otherwise.
  ///
  bool checkDurationSplines() const {
    std::size_t i = 0;
    bool is_positive = true;
    while (is_positive && i < duration_splines_.size()) {
      is_positive = (duration_splines_.at(i) > 0.);
      i++;
    }
    return is_positive;
  }
  /*Operations*/

  /*Helpers*/
 public:
  /// \brief Get dimension of curve.
  /// \return dimension of curve.
  std::size_t virtual dim() const { return dim_; }
  /// \brief Get the minimum time for which the curve is defined
  /// \return \f$t_{min}\f$, lower bound of time range.
  Time virtual min() const { return time_control_points_.front(); }
  /// \brief Get the maximum time for which the curve is defined.
  /// \return \f$t_{max}\f$, upper bound of time range.
  Time virtual max() const { return time_control_points_.back(); }
  /// \brief Get the degree of the curve.
  /// \return \f$degree\f$, the degree of the curve.
  virtual std::size_t degree() const { return degree_; }
  /*Helpers*/

  /*Attributes*/
  /// Dim of curve
  std::size_t dim_;
  /// Vector of pair < Point, Tangent >.
  t_pair_point_tangent_t control_points_;
  /// Vector of Time corresponding to time of each N control points : time at
  /// \f$P_0, P_1, P_2, ..., P_N\f$. Exemple : \f$( 0., 0.5, 0.9, ..., 4.5 )\f$
  /// with values corresponding to times for \f$P_0, P_1, P_2, ..., P_N\f$
  /// respectively.
  vector_time_t time_control_points_;
  /// Vector of Time corresponding to time duration of each subspline.<br>
  /// For N control points with time \f$T_{P_0}, T_{P_1}, T_{P_2}, ...,
  /// T_{P_N}\f$ respectively, duration of each subspline is : (
  /// T_{P_1}-T_{P_0}, T_{P_2}-T_{P_1}, ..., T_{P_N}-T_{P_{N-1} )<br> It
  /// contains \f$N-1\f$ durations.
  vector_time_t duration_splines_;
  /// Starting time of cubic hermite spline : T_min_ is equal to first time of
  /// control points.
  /*const*/ Time T_min_;
  /// Ending time of cubic hermite spline : T_max_ is equal to last time of
  /// control points.
  /*const*/ Time T_max_;
  /// Number of control points (pairs).
  std::size_t size_;
  /// Degree (Cubic so degree 3)
  std::size_t degree_;
  /*Attributes*/

  // Serialization of the class
  friend class boost::serialization::access;

  template <class Archive>
  void serialize(Archive& ar, const unsigned int version) {
    if (version) {
      // Do something depending on version ?
    }
    ar& BOOST_SERIALIZATION_BASE_OBJECT_NVP(curve_abc_t);
    ar& boost::serialization::make_nvp("dim", dim_);
    ar& boost::serialization::make_nvp("control_points", control_points_);
    ar& boost::serialization::make_nvp("time_control_points",
                                       time_control_points_);
    ar& boost::serialization::make_nvp("duration_splines", duration_splines_);
    ar& boost::serialization::make_nvp("T_min", T_min_);
    ar& boost::serialization::make_nvp("T_max", T_max_);
    ar& boost::serialization::make_nvp("size", size_);
    ar& boost::serialization::make_nvp("degree", degree_);
  }
};  // End struct Cubic hermite spline
}  // namespace ndcurves

DEFINE_CLASS_TEMPLATE_VERSION(
    SINGLE_ARG(typename Time, typename Numeric, bool Safe, typename Point),
    SINGLE_ARG(ndcurves::cubic_hermite_spline<Time, Numeric, Safe, Point>))
#endif  //_CLASS_CUBICHERMITESPLINE
