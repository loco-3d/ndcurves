/**
 * \file piecewise_curve.h
 * \brief class allowing to create a piecewise curve.
 * \author Jason C.
 * \date 05/2019
 */

#ifndef _CLASS_PIECEWISE_CURVE
#define _CLASS_PIECEWISE_CURVE

#include <boost/serialization/vector.hpp>
#include <fstream>
#include <memory>
#include <sstream>

#include "curve_abc.h"
#include "curve_conversion.h"

namespace ndcurves {
/// \class PiecewiseCurve.
/// \brief Represent a piecewise curve. We can add some new curve,
///        but the starting time of the curve to add should be equal to the
///        ending time of the actual piecewise_curve.<br>\ Example : A piecewise
///        curve composed of three curves cf0, cf1 and cf2 where cf0 is defined
///        between \f$[T0_{min},T0_{max}]\f$, cf1 between
///        \f$[T0_{max},T1_{max}]\f$ and cf2 between \f$[T1_{max},T2_{max}]\f$.
///        On the piecewise polynomial curve, cf0 is located between
///        \f$[T0_{min},T0_{max}[\f$, cf1 between \f$[T0_{max},T1_{max}[\f$ and
///        cf2 between \f$[T1_{max},T2_{max}]\f$.
///
template <typename Time = double, typename Numeric = Time, bool Safe = false,
          typename Point = Eigen::Matrix<Numeric, Eigen::Dynamic, 1>,
          typename Point_derivate = Point,
          typename CurveType =
              curve_abc<Time, Numeric, Safe, Point, Point_derivate> >
struct piecewise_curve
    : public curve_abc<Time, Numeric, Safe, Point, Point_derivate> {
  typedef Point point_t;
  typedef Point_derivate point_derivate_t;
  typedef std::vector<point_t, Eigen::aligned_allocator<point_t> > t_point_t;
  typedef std::vector<point_derivate_t,
                      Eigen::aligned_allocator<point_derivate_t> >
      t_point_derivate_t;
  typedef Time time_t;
  typedef Numeric num_t;
  typedef curve_abc<Time, Numeric, Safe, point_t, point_derivate_t>
      base_curve_t;           // parent class
  typedef CurveType curve_t;  // contained curves base class
  typedef std::shared_ptr<curve_t> curve_ptr_t;
  typedef typename std::vector<curve_ptr_t> t_curve_ptr_t;
  typedef typename std::vector<Time> t_time_t;
  typedef piecewise_curve<Time, Numeric, Safe, Point, Point_derivate, CurveType>
      piecewise_curve_t;
  typedef piecewise_curve<Time, Numeric, Safe, Point_derivate, Point_derivate,
                          typename CurveType::curve_derivate_t>
      piecewise_curve_derivate_t;
  typedef std::shared_ptr<typename piecewise_curve_derivate_t::curve_t>
      curve_derivate_ptr_t;

 public:
  /// \brief Empty constructor. Add at least one curve to call other class
  /// functions.
  ///
  piecewise_curve() : dim_(0), size_(0), T_min_(0), T_max_(0) {}

  /// \brief Constructor.
  /// Initialize a piecewise curve by giving the first curve.
  /// \param cf   : a curve.
  ///
  piecewise_curve(const curve_ptr_t& cf)
      : dim_(0), size_(0), T_min_(0), T_max_(0) {
    add_curve_ptr(cf);
  }

  piecewise_curve(const t_curve_ptr_t& curves_list)
      : dim_(0), size_(0), T_min_(0), T_max_(0) {
    for (typename t_curve_ptr_t::const_iterator it = curves_list.begin();
         it != curves_list.end(); ++it) {
      add_curve_ptr(*it);
    }
  }

  piecewise_curve(const piecewise_curve& other)
      : dim_(other.dim_),
        curves_(other.curves_),
        time_curves_(other.time_curves_),
        size_(other.size_),
        T_min_(other.T_min_),
        T_max_(other.T_max_) {}

  virtual ~piecewise_curve() {}

  virtual point_t operator()(const Time t) const {
    check_if_not_empty();
    if (Safe & !(T_min_ <= t && t <= T_max_)) {
      // std::cout<<"[Min,Max]=["<<T_min_<<","<<T_max_<<"]"<<"
      // t="<<t<<std::endl;
      throw std::out_of_range("can't evaluate piecewise curve, out of range");
    }
    return (*curves_.at(find_interval(t)))(t);
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
      const piecewise_curve_t& other,
      const Numeric prec = Eigen::NumTraits<Numeric>::dummy_precision()) const {
    if (num_curves() != other.num_curves()) return false;
    for (size_t i = 0; i < num_curves(); ++i) {
      if (!curve_at_index(i)->isApprox(other.curve_at_index(i).get(), prec))
        return false;
    }
    return true;
  }

  virtual bool isApprox(
      const base_curve_t* other,
      const Numeric prec = Eigen::NumTraits<Numeric>::dummy_precision()) const {
    const piecewise_curve_t* other_cast =
        dynamic_cast<const piecewise_curve_t*>(other);
    if (other_cast)
      return isApprox(*other_cast, prec);
    else
      return false;
  }

  virtual bool operator==(const piecewise_curve_t& other) const {
    return isApprox(other);
  }

  virtual bool operator!=(const piecewise_curve_t& other) const {
    return !(*this == other);
  }

  ///  \brief Evaluate the derivative of order N of curve at time t.
  ///  \param t : time when to evaluate the spline.
  ///  \param order : order of derivative.
  ///  \return \f$\frac{d^Np(t)}{dt^N}\f$ point corresponding on derivative
  ///  spline of order N at time t.
  ///
  virtual point_derivate_t derivate(const Time t,
                                    const std::size_t order) const {
    check_if_not_empty();
    if (Safe & !(T_min_ <= t && t <= T_max_)) {
      throw std::invalid_argument(
          "can't evaluate piecewise curve, out of range");
    }
    return (*curves_.at(find_interval(t))).derivate(t, order);
  }

  /**
   * @brief compute_derivate return a piecewise_curve which is the derivative of
   * this at given order
   * @param order order of derivative
   * @return
   */
  piecewise_curve_derivate_t* compute_derivate_ptr(
      const std::size_t order) const {
    piecewise_curve_derivate_t* res(new piecewise_curve_derivate_t());
    for (typename t_curve_ptr_t::const_iterator itc = curves_.begin();
         itc < curves_.end(); ++itc) {
      curve_derivate_ptr_t ptr((*itc)->compute_derivate_ptr(order));
      res->add_curve_ptr(ptr);
    }
    return res;
  }

  template <typename Curve>
  void add_curve(const Curve& curve) {
    curve_ptr_t curve_ptr = std::make_shared<Curve>(curve);
    add_curve_ptr(curve_ptr);
  }

  ///  \brief Add a new curve to piecewise curve, which should be defined in
  ///  \f$[T_{min},T_{max}]\f$ where \f$T_{min}\f$
  ///         is equal to \f$T_{max}\f$ of the actual piecewise curve. The curve
  ///         added should be of type Curve as defined in the template.
  ///  \param cf : curve to add.
  ///
  void add_curve_ptr(const curve_ptr_t& cf) {
    if (size_ == 0) {  // first curve added
      dim_ = cf->dim();
    }
    // Check time continuity : Beginning time of cf must be equal to T_max_ of
    // actual piecewise curve.
    if (size_ != 0 && !(fabs(cf->min() - T_max_) < MARGIN)) {
      std::stringstream ss;
      ss << "Can not add new Polynom to PiecewiseCurve : time discontinuity "
            "between T_max_ and pol.min(). Current "
            "T_max is "
         << T_max_ << " new curve min is " << cf->min();
      throw std::invalid_argument(ss.str().c_str());
    }
    if (cf->dim() != dim_) {
      std::stringstream ss;
      ss << "All the curves in a piecewiseCurve should have the same "
            "dimension. Current dim is "
         << dim_ << " dim of the new curve is " << cf->dim();
      throw std::invalid_argument(ss.str().c_str());
    }
    curves_.push_back(cf);
    size_ = curves_.size();
    T_max_ = cf->max();
    if (size_ == 1) {
      // First curve added
      time_curves_.push_back(cf->min());
      T_min_ = cf->min();
    }
    time_curves_.push_back(T_max_);
  }

  ///  \brief Check if the curve is continuous of order given.
  ///  \param order : order of continuity we want to check.
  ///  \return True if the curve is continuous of order given.
  ///
  bool is_continuous(const std::size_t order) {
    check_if_not_empty();
    bool isContinuous = true;
    std::size_t i = 0;
    if (order == 0) {
      point_t value_end, value_start;
      while (isContinuous && i < (size_ - 1)) {
        curve_ptr_t current = curves_.at(i);
        curve_ptr_t next = curves_.at(i + 1);
        value_end = (*current)(current->max());
        value_start = (*next)(next->min());
        if (!value_end.isApprox(value_start, MARGIN)) {
          isContinuous = false;
        }
        i++;
      }
    } else {
      point_derivate_t value_end, value_start;
      while (isContinuous && i < (size_ - 1)) {
        curve_ptr_t current = curves_.at(i);
        curve_ptr_t next = curves_.at(i + 1);
        value_end = current->derivate(current->max(), order);
        value_start = next->derivate(next->min(), order);
        if (!value_end.isApprox(value_start, MARGIN)) {
          isContinuous = false;
        }
        i++;
      }
    }
    return isContinuous;
  }

  /// \brief Get number of curves in piecewise curve.
  /// \return Number of curves in piecewise curve.
  std::size_t num_curves() const { return curves_.size(); }

  /// \brief Get curve corresponding to time t in piecewise curve.
  /// Example : A piecewise curve PC made of two curves : c1 for t in [0,1] and
  /// c2 for t in ]1,2].
  ///           PC.curve_at_time(0.5) will return c1.
  /// \param t : time to select curve.
  /// \return Curve corresponding to time t in piecewise curve.
  curve_ptr_t curve_at_time(const time_t t) const {
    return curves_[find_interval(t)];
  }

  /// \brief Get curve at specified index in piecewise curve.
  /// \param idx : Index of curve to return, from 0 to num_curves-1.
  /// \return curve corresonding to index in piecewise curve.
  curve_ptr_t curve_at_index(const std::size_t idx) const {
    if (Safe && idx >= num_curves()) {
      throw std::length_error(
          "curve_at_index: requested index greater than number of curves in "
          "piecewise_curve instance");
    }
    return curves_[idx];
  }

  /// \brief Convert all curves in piecewise curve into bezier curves.
  /// \return piecewise bezier curve.
  ///
  template <typename Bezier>
  piecewise_curve_t convert_piecewise_curve_to_bezier() {
    check_if_not_empty();
    // check if given Bezier curve have the correct dimension :
    BOOST_STATIC_ASSERT(
        boost::is_same<typename Bezier::point_t, point_t>::value);
    BOOST_STATIC_ASSERT(boost::is_same<typename Bezier::point_derivate_t,
                                       point_derivate_t>::value);
    // Create piecewise curve
    piecewise_curve_t pc_res;
    // Convert and add all other curves (segments)
    for (std::size_t i = 0; i < size_; i++) {
      pc_res.add_curve(bezier_from_curve<Bezier>(*curves_.at(i)));
    }
    return pc_res;
  }

  /// \brief Convert all curves in piecewise curve into cubic hermite curves.
  /// Curves need to be of degree inferior or equal to three.
  /// \return piecewise cubic hermite curve.
  ///
  template <typename Hermite>
  piecewise_curve_t convert_piecewise_curve_to_cubic_hermite() {
    check_if_not_empty();
    // check if given Hermite curve have the correct dimension :
    BOOST_STATIC_ASSERT(
        boost::is_same<typename Hermite::point_t, point_t>::value);
    BOOST_STATIC_ASSERT(boost::is_same<typename Hermite::point_derivate_t,
                                       point_derivate_t>::value);
    // Create piecewise curve
    piecewise_curve_t pc_res;
    // Convert and add all other curves (segments)
    for (std::size_t i = 0; i < size_; i++) {
      pc_res.add_curve(hermite_from_curve<Hermite>(*curves_.at(i)));
    }
    return pc_res;
  }

  /// \brief Convert all curves in piecewise curve into polynomial curves.
  /// \return piecewise polynomial curve.
  ///
  template <typename Polynomial>
  piecewise_curve_t convert_piecewise_curve_to_polynomial() {
    check_if_not_empty();
    // check if given Polynomial curve have the correct dimension :
    BOOST_STATIC_ASSERT(
        boost::is_same<typename Polynomial::point_t, point_t>::value);
    BOOST_STATIC_ASSERT(boost::is_same<typename Polynomial::point_derivate_t,
                                       point_derivate_t>::value);
    // Create piecewise curve
    piecewise_curve_t pc_res;
    // Convert and add all other curves (segments)
    for (std::size_t i = 0; i < size_; i++) {
      pc_res.add_curve(polynomial_from_curve<Polynomial>(*curves_.at(i)));
    }
    return pc_res;
  }

  /// \brief Convert discrete points into piecewise polynomial curve with C0
  /// continuity. \param points : discrete points to convert. \param time_points
  /// : time corresponding to each point in piecewise curve. \return piecewise
  /// polynomial curve of C0 continuity.
  ///
  template <typename Polynomial>
  static piecewise_curve_t convert_discrete_points_to_polynomial(
      t_point_t points, t_time_t time_points) {
    if (Safe & !(points.size() > 1)) {
      // std::cout<<"[Min,Max]=["<<T_min_<<","<<T_max_<<"]"<<"
      // t="<<t<<std::endl;
      throw std::invalid_argument(
          "piecewise_curve::convert_discrete_points_to_polynomial: Error, less "
          "than 2 discrete points");
    }
    if (points.size() != time_points.size()) {
      throw std::invalid_argument(
          "piecewise_curve::convert_discrete_points_to_polynomial: Error, "
          "points and time_points must have the same "
          "size.");
    }
    // check if given Polynomial curve have the correct dimension :
    BOOST_STATIC_ASSERT(
        boost::is_same<typename Polynomial::point_t, point_t>::value);
    BOOST_STATIC_ASSERT(boost::is_same<typename Polynomial::point_derivate_t,
                                       point_derivate_t>::value);
    piecewise_curve_t piecewise_res;

    for (size_t i = 1; i < points.size(); ++i) {
      piecewise_res.add_curve(Polynomial(points[i - 1], points[i],
                                         time_points[i - 1], time_points[i]));
    }
    return piecewise_res;
  }

  /// \brief Convert discrete points into piecewise polynomial curve with C1
  /// continuity. \param points : discrete points to convert. \param
  /// points_derivative : derivative of order 1 corresponding to each point in
  /// piecewise curve. \param time_points : time corresponding to each point in
  /// piecewise curve. \return piecewise polynomial curve of C1 continuity.
  ///
  template <typename Polynomial>
  static piecewise_curve_t convert_discrete_points_to_polynomial(
      t_point_t points, t_point_derivate_t points_derivative,
      t_time_t time_points) {
    if (Safe & !(points.size() > 1)) {
      // std::cout<<"[Min,Max]=["<<T_min_<<","<<T_max_<<"]"<<"
      // t="<<t<<std::endl;
      throw std::invalid_argument(
          "piecewise_curve::convert_discrete_points_to_polynomial: Error, less "
          "than 2 discrete points");
    }
    if (points.size() != time_points.size()) {
      throw std::invalid_argument(
          "piecewise_curve::convert_discrete_points_to_polynomial: Error, "
          "points and time_points must have the same "
          "size.");
    }
    if (points.size() != points_derivative.size()) {
      throw std::invalid_argument(
          "piecewise_curve::convert_discrete_points_to_polynomial: Error, "
          "points and points_derivative must have the "
          "same size.");
    }
    // check if given Polynomial curve have the correct dimension :
    BOOST_STATIC_ASSERT(
        boost::is_same<typename Polynomial::point_t, point_t>::value);
    BOOST_STATIC_ASSERT(boost::is_same<typename Polynomial::point_derivate_t,
                                       point_derivate_t>::value);
    piecewise_curve_t piecewise_res;

    for (size_t i = 1; i < points.size(); ++i) {
      piecewise_res.add_curve(
          Polynomial(points[i - 1], points_derivative[i - 1], points[i],
                     points_derivative[i], time_points[i - 1], time_points[i]));
    }
    return piecewise_res;
  }

  /// \brief Convert discrete points into piecewise polynomial curve with C2
  /// continuity. \param points : discrete points to convert. \param
  /// points_derivative : derivative of order 1 corresponding to each point in
  /// piecewise curve. \param points_second_derivative : derivative of order 2
  /// corresponding to each point in piecewise curve. \param time_points : time
  /// corresponding to each point in piecewise curve. \return piecewise
  /// polynomial curve of C2 continuity.
  ///
  template <typename Polynomial>
  static piecewise_curve_t convert_discrete_points_to_polynomial(
      t_point_t points, t_point_derivate_t points_derivative,
      t_point_derivate_t points_second_derivative, t_time_t time_points) {
    if (Safe & !(points.size() > 1)) {
      // std::cout<<"[Min,Max]=["<<T_min_<<","<<T_max_<<"]"<<"
      // t="<<t<<std::endl;
      throw std::invalid_argument(
          "piecewise_curve::convert_discrete_points_to_polynomial: Error, less "
          "than 2 discrete points");
    }
    if (points.size() != time_points.size()) {
      throw std::invalid_argument(
          "piecewise_curve::convert_discrete_points_to_polynomial: Error, "
          "points and time_points must have the same "
          "size.");
    }
    if (points.size() != points_derivative.size()) {
      throw std::invalid_argument(
          "piecewise_curve::convert_discrete_points_to_polynomial: Error, "
          "points and points_derivative must have the "
          "same size.");
    }
    if (points.size() != points_second_derivative.size()) {
      throw std::invalid_argument(
          "piecewise_curve::convert_discrete_points_to_polynomial: Error, "
          "points and points_second_derivative must "
          "have the same size.");
    }
    // check if given Polynomial curve have the correct dimension :
    BOOST_STATIC_ASSERT(
        boost::is_same<typename Polynomial::point_t, point_t>::value);
    BOOST_STATIC_ASSERT(boost::is_same<typename Polynomial::point_derivate_t,
                                       point_derivate_t>::value);
    piecewise_curve_t piecewise_res;

    for (size_t i = 1; i < points.size(); ++i) {
      piecewise_res.add_curve(Polynomial(
          points[i - 1], points_derivative[i - 1],
          points_second_derivative[i - 1], points[i], points_derivative[i],
          points_second_derivative[i], time_points[i - 1], time_points[i]));
    }
    return piecewise_res;
  }

  /**
   * @brief load_piecewise_from_text_file build a piecewise polynomial from a
   * list of discrete points read from a file. The file should contains one
   * points per line, optionally with it's derivative and second derivatives.
   * Each lines should then contains dim, 2*dim or 3*dim values
   * @param filename the (absolute) name of the file to load
   * @param dt the time step between each points in the file
   * @param dim the dimension of the curve
   * @return a piecewise curves containing polynomial connectiong all the points
   * in the file
   */
  template <typename Polynomial>
  static piecewise_curve_t load_piecewise_from_text_file(
      const std::string& filename, const time_t dt, const size_t dim) {
    if (dim <= 0)
      throw std::invalid_argument("The dimension should be strictly positive.");
    if (dt <= 0.)
      throw std::invalid_argument("The time step should be strictly positive.");

    piecewise_curve_t piecewise_res;
    std::ifstream file;
    file.open(filename.c_str());
    point_t last_pos = point_t::Zero(dim), last_vel = point_t::Zero(dim),
            last_acc = point_t::Zero(dim), new_pos = point_t::Zero(dim),
            new_vel = point_t::Zero(dim), new_acc = point_t::Zero(dim);
    bool use_vel, use_acc;
    std::string line;
    // read first line to found out if we use velocity / acceleration :
    std::getline(file, line);
    std::istringstream iss_length(line);
    const size_t length =
        std::distance(std::istream_iterator<std::string>(iss_length),
                      std::istream_iterator<std::string>());
    if (length == dim) {
      use_vel = false;
      use_acc = false;
    } else if (length == dim * 2) {
      use_vel = true;
      use_acc = false;
    } else if (length == dim * 3) {
      use_vel = true;
      use_acc = true;
    } else {
      std::stringstream error;
      error << "The first line of the file shold contains either " << dim
            << ", " << dim * 2 << " or " << dim * 3
            << "values, got : " << length;
      throw std::invalid_argument(error.str());
    }
    // initialize the first points of the trajectory:
    num_t val;
    std::istringstream iss(line);
    for (size_t i = 0; i < dim; ++i) {
      iss >> val;
      last_pos[i] = val;
    }
    if (use_vel) {
      for (size_t i = 0; i < dim; ++i) {
        iss >> val;
        last_vel[i] = val;
      }
    }
    if (use_acc) {
      for (size_t i = 0; i < dim; ++i) {
        iss >> val;
        last_acc[i] = val;
      }
    }

    size_t current_length;
    size_t line_id = 0;
    // parse all lines of the file:
    while (std::getline(file, line)) {
      ++line_id;
      std::istringstream iss_length(line);
      current_length =
          std::distance(std::istream_iterator<std::string>(iss_length),
                        std::istream_iterator<std::string>());
      if (current_length != length) {
        std::stringstream error;
        error << "Cannot parse line " << line_id << " got " << current_length
              << " values instead of " << length;
        throw std::invalid_argument(error.str());
      }
      std::istringstream iss(line);
      // parse the points values from the file:
      for (size_t i = 0; i < dim; ++i) {
        iss >> val;
        new_pos[i] = val;
      }
      if (use_vel) {
        for (size_t i = 0; i < dim; ++i) {
          iss >> val;
          new_vel[i] = val;
        }
      }
      if (use_acc) {
        for (size_t i = 0; i < dim; ++i) {
          iss >> val;
          new_acc[i] = val;
        }
      }
      // append a new curves connectiong this points
      if (use_acc) {
        piecewise_res.add_curve(
            Polynomial(last_pos, last_vel, last_acc, new_pos, new_vel, new_acc,
                       dt * static_cast<time_t>(line_id - 1),
                       dt * static_cast<time_t>(line_id)));
      } else if (use_vel) {
        piecewise_res.add_curve(
            Polynomial(last_pos, last_vel, new_pos, new_vel,
                       dt * static_cast<time_t>(line_id - 1),
                       dt * static_cast<time_t>(line_id)));
      } else {
        piecewise_res.add_curve(
            Polynomial(last_pos, new_pos, dt * static_cast<time_t>(line_id - 1),
                       dt * static_cast<time_t>(line_id)));
      }
      last_pos = new_pos;
      last_vel = new_vel;
      last_acc = new_acc;
    }

    file.close();
    return piecewise_res;
  }

 private:
  /// \brief Get index of the interval corresponding to time t for the
  /// interpolation. \param t : time where to look for interval. \return Index
  /// of interval for time t.
  ///
  std::size_t find_interval(const Numeric t) const {
    // time before first control point time.
    if (t < time_curves_[0]) {
      return 0;
    }
    // time is after last control point time
    if (t > time_curves_[size_ - 1]) {
      return size_ - 1;
    }

    std::size_t left_id = 0;
    std::size_t right_id = size_ - 1;
    while (left_id <= right_id) {
      const std::size_t middle_id = left_id + (right_id - left_id) / 2;
      if (time_curves_.at(middle_id) < t) {
        left_id = middle_id + 1;
      } else if (time_curves_.at(middle_id) > t) {
        right_id = middle_id - 1;
      } else {
        return middle_id;
      }
    }
    return left_id - 1;
  }

  void check_if_not_empty() const {
    if (curves_.size() == 0) {
      throw std::runtime_error("Error in piecewise curve : No curve added");
    }
  }

  /*Helpers*/
 public:
  /// \brief Get dimension of curve.
  /// \return dimension of curve.
  std::size_t virtual dim() const { return dim_; };
  /// \brief Get the minimum time for which the curve is defined
  /// \return \f$t_{min}\f$, lower bound of time range.
  Time virtual min() const { return T_min_; }
  /// \brief Get the maximum time for which the curve is defined.
  /// \return \f$t_{max}\f$, upper bound of time range.
  Time virtual max() const { return T_max_; }
  /// \brief Get the degree of the curve.
  /// \return \f$degree\f$, the degree of the curve.
  virtual std::size_t degree() const {
    throw std::runtime_error(
        "degree() method is not implemented for this type of curve.");
  }
  std::size_t getNumberCurves() { return curves_.size(); }
  /*Helpers*/

  /* Attributes */
  std::size_t dim_;       // Dim of curve
  t_curve_ptr_t curves_;  // for curves 0/1/2 : [ curve0, curve1, curve2 ]
  t_time_t time_curves_;  // for curves 0/1/2 : [ Tmin0, Tmax0,Tmax1,Tmax2 ]
  std::size_t size_;  // Number of segments in piecewise curve = size of curves_
  Time T_min_, T_max_;
  /* Attributes */

  // Serialization of the class
  friend class boost::serialization::access;

  template <class Archive>
  void serialize(Archive& ar, const unsigned int version) {
    if (version) {
      // Do something depending on version ?
    }
    ar& BOOST_SERIALIZATION_BASE_OBJECT_NVP(base_curve_t);
    ar& boost::serialization::make_nvp("dim", dim_);
    ar& boost::serialization::make_nvp("curves", curves_);
    ar& boost::serialization::make_nvp("time_curves", time_curves_);
    ar& boost::serialization::make_nvp("size", size_);
    ar& boost::serialization::make_nvp("T_min", T_min_);
    ar& boost::serialization::make_nvp("T_max", T_max_);
  }
};  // End struct piecewise curve
}  // namespace ndcurves

DEFINE_CLASS_TEMPLATE_VERSION(
    SINGLE_ARG(typename Time, typename Numeric, bool Safe, typename Point,
               typename Point_derivate, typename CurveType),
    SINGLE_ARG(ndcurves::piecewise_curve<Time, Numeric, Safe, Point,
                                         Point_derivate, CurveType>))

#endif  // _CLASS_PIECEWISE_CURVE
