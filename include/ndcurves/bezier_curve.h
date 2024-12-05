/**
 * \file bezier_curve.h
 * \brief class allowing to create a Bezier curve of dimension 1 <= n <= 3.
 * \author Steve T.
 * \version 0.1
 * \date 06/17/2013
 */

#ifndef _CLASS_BEZIERCURVE
#define _CLASS_BEZIERCURVE

#include <iostream>
#include <stdexcept>
#include <vector>

#include "MathDefs.h"
#include "bernstein.h"
#include "cross_implementation.h"
#include "curve_abc.h"
#include "curve_constraint.h"
#include "piecewise_curve.h"

namespace ndcurves {
/// \class BezierCurve.
/// \brief Represents a Bezier curve of arbitrary dimension and order.
/// For degree lesser than 4, the evaluation is analitycal. Otherwise
/// the bernstein polynoms are used to evaluate the spline at a given location.
///
template <typename Time = double, typename Numeric = Time, bool Safe = false,
          typename Point = Eigen::Matrix<Numeric, Eigen::Dynamic, 1> >
struct bezier_curve : public curve_abc<Time, Numeric, Safe, Point> {
  typedef Point point_t;
  typedef Eigen::Matrix<Numeric, Eigen::Dynamic, 1> vector_x_t;
  typedef Eigen::Ref<const vector_x_t> vector_x_ref_t;
  typedef Time time_t;
  typedef Numeric num_t;
  typedef curve_constraints<point_t> curve_constraints_t;
  typedef std::vector<point_t, Eigen::aligned_allocator<point_t> > t_point_t;
  typedef typename t_point_t::const_iterator cit_point_t;
  typedef bezier_curve<Time, Numeric, Safe, Point> bezier_curve_t;
  typedef std::shared_ptr<bezier_curve_t> bezier_curve_ptr_t;
  typedef piecewise_curve<Time, Numeric, Safe, point_t, point_t, bezier_curve_t>
      piecewise_curve_t;
  typedef curve_abc<Time, Numeric, Safe, point_t> curve_abc_t;  // parent class
  typedef typename curve_abc_t::curve_ptr_t curve_ptr_t;

  /* Constructors - destructors */
 public:
  /// \brief Empty constructor. Curve obtained this way can not perform other
  /// class functions.
  ///
  bezier_curve() : dim_(0), T_min_(0), T_max_(0) {}

  /// \brief Constructor.
  /// Given the first and last point of a control points set, create the bezier
  /// curve. \param PointsBegin   : an iterator pointing to the first element of
  /// a control point container. \param PointsEnd     : an iterator pointing to
  /// the last element of a control point container. \param T_min         :
  /// lower bound of time, curve will be defined for time in [T_min, T_max].
  /// \param T_max         : upper bound of time, curve will be defined for time
  /// in [T_min, T_max]. \param mult_T        : ... (default value is 1.0).
  ///
  template <typename In>
  bezier_curve(In PointsBegin, In PointsEnd, const time_t T_min = 0.,
               const time_t T_max = 1., const time_t mult_T = 1.)
      : dim_(PointsBegin->size()),
        T_min_(T_min),
        T_max_(T_max),
        mult_T_(mult_T),
        size_(std::distance(PointsBegin, PointsEnd)),
        degree_(size_ - 1),
        bernstein_(ndcurves::makeBernstein<num_t>((unsigned int)degree_)) {
    if (bernstein_.size() != size_) {
      throw std::invalid_argument("Invalid size of polynomial");
    }
    In it(PointsBegin);
    if (Safe && (size_ < 1 || T_max_ <= T_min_)) {
      throw std::invalid_argument(
          "can't create bezier min bound is higher than max bound");
    }
    for (; it != PointsEnd; ++it) {
      if (Safe && static_cast<size_t>(it->size()) != dim_)
        throw std::invalid_argument(
            "All the control points must have the same dimension.");
      control_points_.push_back(*it);
    }
  }

  /// \brief Constructor with constraints.
  /// This constructor will add 4 points (2 after the first one, 2 before the
  /// last one) to ensure that velocity and acceleration constraints are
  /// respected. \param PointsBegin   : an iterator pointing to the first
  /// element of a control point container. \param PointsEnd     : an iterator
  /// pointing to the last element of a control point container. \param
  /// constraints : constraints applying on start / end velocities and
  /// acceleration. \param T_min         : lower bound of time, curve will be
  /// defined for time in [T_min, T_max]. \param T_max         : upper bound of
  /// time, curve will be defined for time in [T_min, T_max]. \param mult_T :
  /// ... (default value is 1.0).
  ///
  template <typename In>
  bezier_curve(In PointsBegin, In PointsEnd,
               const curve_constraints_t& constraints, const time_t T_min = 0.,
               const time_t T_max = 1., const time_t mult_T = 1.)
      : dim_(PointsBegin->size()),
        T_min_(T_min),
        T_max_(T_max),
        mult_T_(mult_T),
        size_(std::distance(PointsBegin, PointsEnd) + 4),
        degree_(size_ - 1),
        bernstein_(ndcurves::makeBernstein<num_t>((unsigned int)degree_)) {
    if (Safe && (size_ < 1 || T_max_ <= T_min_)) {
      throw std::invalid_argument(
          "can't create bezier min bound is higher than max bound");
    }
    t_point_t updatedList =
        add_constraints<In>(PointsBegin, PointsEnd, constraints);
    for (cit_point_t cit = updatedList.begin(); cit != updatedList.end();
         ++cit) {
      if (Safe && static_cast<size_t>(cit->size()) != dim_)
        throw std::invalid_argument(
            "All the control points must have the same dimension.");
      control_points_.push_back(*cit);
    }
  }

  bezier_curve(const bezier_curve& other)
      : dim_(other.dim_),
        T_min_(other.T_min_),
        T_max_(other.T_max_),
        mult_T_(other.mult_T_),
        size_(other.size_),
        degree_(other.degree_),
        bernstein_(other.bernstein_),
        control_points_(other.control_points_) {}

  ///\brief Destructor
  virtual ~bezier_curve() {}

  /*Operations*/
  ///  \brief Evaluation of the bezier curve at time t.
  ///  \param t : time when to evaluate the curve.
  ///  \return \f$x(t)\f$ point corresponding on curve at time t.
  virtual point_t operator()(const time_t t) const {
    check_conditions();
    if (Safe & !(T_min_ <= t && t <= T_max_)) {
      throw std::invalid_argument(
          "can't evaluate bezier curve, time t is out of range");  // TODO
    }
    if (size_ == 1) {
      return mult_T_ * control_points_[0];
    } else {
      return evalHorner(t);
    }
  }

  /**
   * @brief isApprox check if other and *this are approximately equals.
   * Only two curves of the same class can be approximately equals, for
   * comparison between different type of curves see isEquivalent
   * @param other the other curve to check
   * @param prec the precision threshold, default
   * Eigen::NumTraits<Numeric>::dummy_precision()
   * @return true if the two curves are approximately equals
   */
  bool isApprox(
      const bezier_curve_t& other,
      const Numeric prec = Eigen::NumTraits<Numeric>::dummy_precision()) const {
    bool equal = ndcurves::isApprox<num_t>(T_min_, other.min()) &&
                 ndcurves::isApprox<num_t>(T_max_, other.max()) &&
                 dim_ == other.dim() && degree_ == other.degree() &&
                 size_ == other.size_ &&
                 ndcurves::isApprox<Numeric>(mult_T_, other.mult_T_) &&
                 bernstein_ == other.bernstein_;
    if (!equal) return false;
    for (size_t i = 0; i < size_; ++i) {
      if (!control_points_.at(i).isApprox(other.control_points_.at(i), prec))
        return false;
    }
    return true;
  }

  virtual bool isApprox(
      const curve_abc_t* other,
      const Numeric prec = Eigen::NumTraits<Numeric>::dummy_precision()) const {
    const bezier_curve_t* other_cast =
        dynamic_cast<const bezier_curve_t*>(other);
    if (other_cast)
      return isApprox(*other_cast, prec);
    else
      return false;
  }

  virtual bool operator==(const bezier_curve_t& other) const {
    return isApprox(other);
  }

  virtual bool operator!=(const bezier_curve_t& other) const {
    return !(*this == other);
  }

  ///  \brief Compute the derived curve at order N.
  ///  Computes the derivative order N, \f$\frac{d^Nx(t)}{dt^N}\f$ of bezier
  ///  curve of parametric equation x(t). \param order : order of derivative.
  ///  \return \f$\frac{d^Nx(t)}{dt^N}\f$ derivative order N of the curve.
  bezier_curve_t compute_derivate(const std::size_t order) const {
    check_conditions();
    if (order == 0) {
      return *this;
    }
    t_point_t derived_wp;
    for (typename t_point_t::const_iterator pit = control_points_.begin();
         pit != control_points_.end() - 1; ++pit) {
      derived_wp.push_back((num_t)degree_ * (*(pit + 1) - (*pit)));
    }
    if (derived_wp.empty()) {
      derived_wp.push_back(point_t::Zero(dim_));
    }
    bezier_curve_t deriv(derived_wp.begin(), derived_wp.end(), T_min_, T_max_,
                         mult_T_ * (1. / (T_max_ - T_min_)));
    return deriv.compute_derivate(order - 1);
  }

  ///  \brief Compute the derived curve at order N.
  ///  \param order : order of derivative.
  ///  \return A pointer to \f$\frac{d^Nx(t)}{dt^N}\f$ derivative order N of the
  ///  curve.
  bezier_curve_t* compute_derivate_ptr(const std::size_t order) const {
    return new bezier_curve_t(compute_derivate(order));
  }

  ///  \brief Compute the primitive of the curve at order N.
  ///  Computes the primitive at order N of bezier curve of parametric equation
  ///  \f$x(t)\f$. <br> At order \f$N=1\f$, the primitve \f$X(t)\f$ of
  ///  \f$x(t)\f$ is such as \f$\frac{dX(t)}{dt} = x(t)\f$. \param order : order
  ///  of the primitive. \param init  : constant valuefor the first point of the
  ///  primitive (can tipycally be zero) \return primitive at order N of x(t).
  bezier_curve_t compute_primitive(const std::size_t order,
                                   const point_t& init) const {
    check_conditions();
    if (order == 0) {
      return *this;
    }
    num_t new_degree_inv = 1. / ((num_t)(degree_ + 1));
    t_point_t n_wp;
    point_t current_sum(init);
    n_wp.push_back(current_sum);
    for (typename t_point_t::const_iterator pit = control_points_.begin();
         pit != control_points_.end(); ++pit) {
      current_sum += *pit;
      n_wp.push_back(current_sum * new_degree_inv);
    }
    bezier_curve_t integ(n_wp.begin(), n_wp.end(), T_min_, T_max_,
                         mult_T_ * (T_max_ - T_min_));
    return integ.compute_primitive(order - 1);
  }

  bezier_curve_t compute_primitive(const std::size_t order) const {
    return compute_primitive(order, point_t::Zero(dim_));
  }

  bezier_curve_t* compute_primitive_ptr(const std::size_t order,
                                        const point_t& init) const {
    return new bezier_curve_t(compute_primitive(order, init));
  }

  ///  \brief Computes a Bezier curve of order degrees higher than the current
  ///  curve, but strictly equivalent. Order elevation is required for addition
  ///  / substraction and other comparison operations. \param order : number of
  ///  order the curve must be updated \return An equivalent Bezier, with one
  ///  more degree.
  bezier_curve_t elevate(const std::size_t order) const {
    t_point_t new_waypoints = control_points_, temp_waypoints;
    for (std::size_t i = 1; i <= order; ++i) {
      num_t new_degree_inv = 1. / ((num_t)(degree_ + i));
      temp_waypoints.push_back(*new_waypoints.begin());
      num_t idx_deg_inv = 0.;
      for (typename t_point_t::const_iterator pit = new_waypoints.begin() + 1;
           pit != new_waypoints.end(); ++pit) {
        idx_deg_inv += new_degree_inv;
        temp_waypoints.push_back(idx_deg_inv * (*(pit - 1)) +
                                 (1 - idx_deg_inv) * (*pit));
      }
      temp_waypoints.push_back(*(new_waypoints.end() - 1));
      new_waypoints = temp_waypoints;
      temp_waypoints.clear();
    }
    return bezier_curve_t(new_waypoints.begin(), new_waypoints.end(), T_min_,
                          T_max_, mult_T_);
  }

  ///  \brief Elevate the Bezier curve of order degrees higher than the current
  ///  curve, but strictly equivalent. Order elevation is required for addition
  ///  / substraction and other comparison operations. \param order : number of
  ///  order the curve must be updated
  void elevate_self(const std::size_t order) {
    if (order > 0) (*this) = elevate(order);
  }

  ///  \brief Evaluate the derivative order N of curve at time t.
  ///  If derivative is to be evaluated several times, it is
  ///  rather recommended to compute derived curve using compute_derivate.
  ///  \param order : order of derivative.
  ///  \param t : time when to evaluate the curve.
  ///  \return \f$\frac{d^Nx(t)}{dt^N}\f$ point corresponding on derived curve
  ///  of order N at time t.
  ///
  virtual point_t derivate(const time_t t, const std::size_t order) const {
    return compute_derivate(order)(t);
  }

  /// \brief Evaluate all Bernstein polynomes for a certain degree.
  /// A bezier curve with N control points is represented by : \f$x(t) =
  /// \sum_{i=0}^{N} B_i^N(t) P_i\f$ with \f$ B_i^N(t) = \binom{N}{i}t^i
  /// (1-t)^{N-i} \f$.<br/> Warning: the horner scheme is about 100 times faster
  /// than this method.<br> This method will probably be removed in the future
  /// as the computation of bernstein polynomial is very costly. \param t : time
  /// when to evaluate the curve. \return \f$x(t)\f$ point corresponding on
  /// curve at time t.
  ///
  point_t evalBernstein(const Numeric t) const {
    const Numeric u = (t - T_min_) / (T_max_ - T_min_);
    point_t res = point_t::Zero(dim_);
    typename t_point_t::const_iterator control_points_it =
        control_points_.begin();
    for (typename std::vector<Bern<Numeric> >::const_iterator cit =
             bernstein_.begin();
         cit != bernstein_.end(); ++cit, ++control_points_it) {
      res += cit->operator()(u) * (*control_points_it);
    }
    return res * mult_T_;
  }

  /// \brief Evaluate all Bernstein polynomes for a certain degree using
  /// Horner's scheme. A bezier curve with N control points is expressed as :
  /// \f$x(t) = \sum_{i=0}^{N} B_i^N(t) P_i\f$.<br> To evaluate the position on
  /// curve at time t,we can apply the Horner's scheme : <br> \f$ x(t) =
  /// (1-t)^N(\sum_{i=0}^{N} \binom{N}{i} \frac{1-t}{t}^i P_i) \f$.<br> Horner's
  /// scheme : for a polynom of degree N expressed by : <br> \f$x(t) = a_0 +
  /// a_1t + a_2t^2 + ... + a_nt^n\f$ where \f$number of additions = N\f$ /
  /// f$number of multiplication = N!\f$<br> Using Horner's method, the polynom
  /// is transformed into : <br> \f$x(t) = a_0 + t(a_1 + t(a_2+t(...))\f$ with N
  /// additions and multiplications. \param t : time when to evaluate the curve.
  /// \return \f$x(t)\f$ point corresponding on curve at time t.
  ///
  point_t evalHorner(const Numeric t) const {
    const Numeric u = (t - T_min_) / (T_max_ - T_min_);
    typename t_point_t::const_iterator control_points_it =
        control_points_.begin();
    Numeric u_op, bc, tn;
    u_op = 1.0 - u;
    bc = 1;
    tn = 1;
    point_t tmp = (*control_points_it) * u_op;
    ++control_points_it;
    for (unsigned int i = 1; i < degree_; i++, ++control_points_it) {
      tn = tn * u;
      bc = bc * ((num_t)(degree_ - i + 1)) / i;
      tmp = (tmp + tn * bc * (*control_points_it)) * u_op;
    }
    return (tmp + tn * u * (*control_points_it)) * mult_T_;
  }

  const t_point_t& waypoints() const { return control_points_; }

  const point_t waypointAtIndex(const std::size_t index) const {
    point_t waypoint;
    if (index < control_points_.size()) {
      waypoint = control_points_[index];
    }
    return waypoint;
  }

  /// \brief Evaluate the curve value at time t using deCasteljau algorithm.
  /// The algorithm will compute the \f$N-1\f$ centroids of parameters
  /// \f${t,1-t}\f$ of consecutive \f$N\f$ control points of bezier curve, and
  /// perform it iteratively until getting one point in the list which will be
  /// the evaluation of bezier curve at time \f$t\f$. \param t : time when to
  /// evaluate the curve. \return \f$x(t)\f$ point corresponding on curve at
  /// time t.
  ///
  point_t evalDeCasteljau(const Numeric t) const {
    // normalize time :
    const Numeric u = (t - T_min_) / (T_max_ - T_min_);
    t_point_t pts = deCasteljauReduction(waypoints(), u);
    while (pts.size() > 1) {
      pts = deCasteljauReduction(pts, u);
    }
    return pts[0] * mult_T_;
  }

  t_point_t deCasteljauReduction(const Numeric t) const {
    const Numeric u = (t - T_min_) / (T_max_ - T_min_);
    return deCasteljauReduction(waypoints(), u);
  }

  /// \brief Compute de Casteljau's reduction of the given list of points at
  /// time t. For the list \f$pts\f$ of N points, compute a new list of points
  /// of size N-1 :<br> \f$<br>( pts[0]*(1-t)+pts[1], pts[1]*(1-t)+pts[2], ...,
  /// pts[0]*(N-2)+pts[N-1] )\f$<br> with t the time when to evaluate bezier
  /// curve.<br>\ The new list contains centroid of parameters \f${t,1-t}\f$ of
  /// consecutive points in the list. \param pts : list of points. \param u   :
  /// NORMALIZED time when to evaluate the curve. \return reduced list of point
  /// (size of pts - 1).
  ///
  t_point_t deCasteljauReduction(const t_point_t& pts, const Numeric u) const {
    if (u < 0 || u > 1) {
      throw std::out_of_range("In deCasteljau reduction : u is not in [0;1]");
    }
    if (pts.size() == 1) {
      return pts;
    }

    t_point_t new_pts;
    for (cit_point_t cit = pts.begin(); cit != (pts.end() - 1); ++cit) {
      new_pts.push_back((1 - u) * (*cit) + u * (*(cit + 1)));
    }
    return new_pts;
  }

  /// \brief Split the bezier curve in 2 at time t.
  /// \param t : list of points.
  /// \return pair containing the first element of both bezier curve obtained.
  ///
  std::pair<bezier_curve_t, bezier_curve_t> split(const Numeric t) const {
    check_conditions();
    if (fabs(t - T_max_) < MARGIN) {
      throw std::runtime_error(
          "can't split curve, interval range is equal to original curve");
    }
    t_point_t wps_first(size_), wps_second(size_);
    const Numeric u = (t - T_min_) / (T_max_ - T_min_);
    t_point_t casteljau_pts = waypoints();
    wps_first[0] = casteljau_pts.front();
    wps_second[degree_] = casteljau_pts.back();
    size_t id = 1;
    while (casteljau_pts.size() > 1) {
      casteljau_pts = deCasteljauReduction(casteljau_pts, u);
      wps_first[id] = casteljau_pts.front();
      wps_second[degree_ - id] = casteljau_pts.back();
      ++id;
    }
    bezier_curve_t c_first(wps_first.begin(), wps_first.end(), T_min_, t,
                           mult_T_);
    bezier_curve_t c_second(wps_second.begin(), wps_second.end(), t, T_max_,
                            mult_T_);
    return std::make_pair(c_first, c_second);
  }

  /// \brief Split the bezier curve in several curves, all accessible
  /// within a piecewise_curve_t.
  /// \param times : list of times of size n.
  /// \return a piecewise_curve_t comprising n+1 curves
  ///
  piecewise_curve_t split(const vector_x_t& times) const {
    std::vector<bezier_curve_t> curves;
    bezier_curve_t current = *this;
    for (int i = 0; i < times.rows(); ++i) {
      std::pair<bezier_curve_t, bezier_curve_t> pairsplit =
          current.split(times[i]);
      curves.push_back(pairsplit.first);
      current = pairsplit.second;
    }
    curves.push_back(current);
    piecewise_curve_t res;
    for (typename std::vector<bezier_curve_t>::const_iterator cit =
             curves.begin();
         cit != curves.end(); ++cit) {
      typename piecewise_curve_t::curve_ptr_t ptr(new bezier_curve_t(*cit));
      res.add_curve_ptr(ptr);
    }
    return res;
  }

  /// \brief Extract a bezier curve defined between \f$[t_1,t_2]\f$ from the
  /// actual bezier curve
  ///        defined between \f$[T_{min},T_{max}]\f$ with \f$T_{min} \leq t_1
  ///        \leq t_2 \leq T_{max}\f$.
  /// \param t1 : start time of bezier curve extracted.
  /// \param t2 : end time of bezier curve extracted.
  /// \return bezier curve extract defined between \f$[t_1,t_2]\f$.
  ///
  bezier_curve_t extract(const Numeric t1, const Numeric t2) {
    if (t1 < T_min_ || t1 > T_max_ || t2 < T_min_ || t2 > T_max_) {
      throw std::out_of_range("In Extract curve : times out of bounds");
    }
    if (fabs(t1 - T_min_) < MARGIN &&
        fabs(t2 - T_max_) < MARGIN)  // t1=T_min and t2=T_max
    {
      return bezier_curve_t(waypoints().begin(), waypoints().end(), T_min_,
                            T_max_, mult_T_);
    }
    if (fabs(t1 - T_min_) < MARGIN)  // t1=T_min
    {
      return split(t2).first;
    }
    if (fabs(t2 - T_max_) < MARGIN)  // t2=T_max
    {
      return split(t1).second;
    }
    std::pair<bezier_curve_t, bezier_curve_t> c_split = this->split(t1);
    return c_split.second.split(t2).first;
  }

  ///  \brief Compute the cross product of the current bezier curve by another
  ///  bezier curve.
  /// The cross product p1Xp2 of 2 bezier curves p1 and p2 is defined such that
  /// forall t, p1Xp2(t) = p1(t) X p2(t), with X designing the cross product.
  /// This method of course only makes sense for dimension 3 curves.
  /// It assumes that a method point_t cross(const point_t&, const point_t&) has
  /// been defined
  ///  \param pOther other polynomial to compute the cross product with.
  ///  \return a new polynomial defining the cross product between this and
  ///  pother
  bezier_curve_t cross(const bezier_curve_t& g) const {
    // See Farouki and Rajan 1988 Alogirthms for polynomials in Bernstein form
    // and http://web.mit.edu/hyperbook/Patrikalakis-Maekawa-Cho/node10.html
    assert_operator_compatible(g);
    if (dim() != 3)
      throw std::invalid_argument(
          "Can't perform cross product on Bezier curves with dimensions != 3 ");
    int m = (int)(degree());
    int n = (int)(g.degree());
    unsigned int mj, n_ij, mn_i;
    t_point_t new_waypoints;
    for (int i = 0; i <= m + n; ++i) {
      bezier_curve_t::point_t current_point =
          bezier_curve_t::point_t::Zero(dim());
      for (int j = std::max(0, i - n); j <= std::min(m, i); ++j) {
        mj = bin(m, j);
        n_ij = bin(n, i - j);
        mn_i = bin(m + n, i);
        num_t mul = num_t(mj * n_ij) / num_t(mn_i);
        current_point +=
            mul * ndcurves::cross(waypointAtIndex(j), g.waypointAtIndex(i - j));
      }
      new_waypoints.push_back(current_point);
    }
    return bezier_curve_t(new_waypoints.begin(), new_waypoints.end(), min(),
                          max(), mult_T_ * g.mult_T_);
  }

  ///  \brief Compute the cross product of the current bezier b by a point
  ///  point.
  /// The cross product pXpoint of is defined such that
  /// forall t, bXpoint(t) = b(t) X point, with X designing the cross product.
  /// This method of course only makes sense for dimension 3 polynomials.
  ///  \param point point to compute the cross product with.
  ///  \return a new polynomial defining the cross product between this and
  ///  point
  bezier_curve_t cross(const bezier_curve_t::point_t& point) const {
    // See Farouki and Rajan 1988 Alogirthms for polynomials in Bernstein form
    // and http://web.mit.edu/hyperbook/Patrikalakis-Maekawa-Cho/node10.html
    if (dim() != 3)
      throw std::invalid_argument(
          "Can't perform cross product on Bezier curves with dimensions != 3 ");
    t_point_t new_waypoints;
    for (typename t_point_t::const_iterator cit = waypoints().begin();
         cit != waypoints().end(); ++cit) {
      new_waypoints.push_back(ndcurves::cross(*cit, point));
    }
    return bezier_curve_t(new_waypoints.begin(), new_waypoints.end(), min(),
                          max(), mult_T_);
  }

  bezier_curve_t& operator+=(const bezier_curve_t& other) {
    assert_operator_compatible(other);
    bezier_curve_t other_elevated =
        other *
        (other.mult_T_ / this->mult_T_);  // TODO remove mult_T_ from Bezier
    if (other.degree() > degree()) {
      elevate_self(other.degree() - degree());
    } else if (other_elevated.degree() < degree()) {
      other_elevated.elevate_self(degree() - other_elevated.degree());
    }
    typename t_point_t::const_iterator otherit =
        other_elevated.control_points_.begin();
    for (typename t_point_t::iterator it = control_points_.begin();
         it != control_points_.end(); ++it, ++otherit) {
      (*it) += (*otherit);
    }
    return *this;
  }

  bezier_curve_t& operator-=(const bezier_curve_t& other) {
    assert_operator_compatible(other);
    bezier_curve_t other_elevated = other * (other.mult_T_ / this->mult_T_);
    if (other.degree() > degree()) {
      elevate_self(other.degree() - degree());
    } else if (other_elevated.degree() < degree()) {
      other_elevated.elevate_self(degree() - other_elevated.degree());
    }
    typename t_point_t::const_iterator otherit =
        other_elevated.control_points_.begin();
    for (typename t_point_t::iterator it = control_points_.begin();
         it != control_points_.end(); ++it, ++otherit) {
      (*it) -= (*otherit);
    }
    return *this;
  }

  bezier_curve_t& operator+=(const bezier_curve_t::point_t& point) {
    for (typename t_point_t::iterator it = control_points_.begin();
         it != control_points_.end(); ++it) {
      (*it) += point;
    }
    return *this;
  }

  bezier_curve_t& operator-=(const bezier_curve_t::point_t& point) {
    for (typename t_point_t::iterator it = control_points_.begin();
         it != control_points_.end(); ++it) {
      (*it) -= point;
    }
    return *this;
  }

  bezier_curve_t& operator/=(const double d) {
    for (typename t_point_t::iterator it = control_points_.begin();
         it != control_points_.end(); ++it) {
      (*it) /= d;
    }
    return *this;
  }

  bezier_curve_t& operator*=(const double d) {
    for (typename t_point_t::iterator it = control_points_.begin();
         it != control_points_.end(); ++it) {
      (*it) *= d;
    }
    return *this;
  }

 private:
  /// \brief Ensure constraints of bezier curve.
  /// Add 4 points (2 after the first one, 2 before the last one) to Bezier
  /// curve to ensure that velocity and acceleration constraints are respected.
  ///
  template <typename In>
  t_point_t add_constraints(In PointsBegin, In PointsEnd,
                            const curve_constraints_t& constraints) {
    t_point_t res;
    num_t T = T_max_ - T_min_;
    num_t T_square = T * T;
    point_t P0, P1, P2, P_n_2, P_n_1, PN;
    P0 = *PointsBegin;
    PN = *(PointsEnd - 1);
    P1 = P0 + constraints.init_vel * T / (num_t)degree_;
    P_n_1 = PN - constraints.end_vel * T / (num_t)degree_;
    P2 = constraints.init_acc * T_square / (num_t)(degree_ * (degree_ - 1)) +
         2 * P1 - P0;
    P_n_2 = constraints.end_acc * T_square / (num_t)(degree_ * (degree_ - 1)) +
            2 * P_n_1 - PN;
    res.push_back(P0);
    res.push_back(P1);
    res.push_back(P2);
    for (In it = PointsBegin + 1; it != PointsEnd - 1; ++it) {
      res.push_back(*it);
    }
    res.push_back(P_n_2);
    res.push_back(P_n_1);
    res.push_back(PN);
    return res;
  }

  void check_conditions() const {
    if (control_points_.size() == 0) {
      throw std::runtime_error(
          "Error in bezier curve : there is no control points set / did you "
          "use empty constructor ?");
    } else if (dim_ == 0) {
      throw std::runtime_error(
          "Error in bezier curve : Dimension of points is zero / did you use "
          "empty constructor ?");
    }
  }

  void assert_operator_compatible(const bezier_curve_t& other) const {
    if ((fabs(min() - other.min()) > MARGIN) ||
        (fabs(max() - other.max()) > MARGIN)) {
      throw std::invalid_argument(
          "Can't perform base operation (+ - ) on two Bezier curves with "
          "different time ranges");
    }
  }

  /*Operations*/

 public:
  /*Helpers*/
  /// \brief Get dimension of curve.
  /// \return dimension of curve.
  std::size_t virtual dim() const { return dim_; };
  /// \brief Get the minimum time for which the curve is defined
  /// \return \f$t_{min}\f$, lower bound of time range.
  virtual time_t min() const { return T_min_; }
  /// \brief Get the maximum time for which the curve is defined.
  /// \return \f$t_{max}\f$, upper bound of time range.
  virtual time_t max() const { return T_max_; }
  /// \brief Get the degree of the curve.
  /// \return \f$degree\f$, the degree of the curve.
  virtual std::size_t degree() const { return degree_; }
  /*Helpers*/

  /* Attributes */
  /// Dim of curve
  std::size_t dim_;
  /// Starting time of cubic hermite spline : T_min_ is equal to first time of
  /// control points.
  /*const*/ time_t T_min_;
  /// Ending time of cubic hermite spline : T_max_ is equal to last time of
  /// control points.
  /*const*/ time_t T_max_;
  /*const*/ time_t mult_T_;
  /*const*/ std::size_t size_;
  /*const*/ std::size_t degree_;
  /*const*/ std::vector<Bern<Numeric> > bernstein_;
  /*const*/ t_point_t control_points_;
  /* Attributes */

  static bezier_curve_t zero(const std::size_t dim, const time_t T = 1.) {
    std::vector<point_t> ts;
    ts.push_back(point_t::Zero(dim));
    return bezier_curve_t(ts.begin(), ts.end(), 0., T);
  }

  // Serialization of the class
  friend class boost::serialization::access;

  template <class Archive>
  void serialize(Archive& ar, const unsigned int version) {
    if (version) {
      // Do something depending on version ?
    }
    ar& BOOST_SERIALIZATION_BASE_OBJECT_NVP(curve_abc_t);
    ar& boost::serialization::make_nvp("dim", dim_);
    ar& boost::serialization::make_nvp("T_min", T_min_);
    ar& boost::serialization::make_nvp("T_max", T_max_);
    ar& boost::serialization::make_nvp("mult_T", mult_T_);
    ar& boost::serialization::make_nvp("size", size_);
    ar& boost::serialization::make_nvp("degree", degree_);
    ar& boost::serialization::make_nvp("bernstein", bernstein_);
    ar& boost::serialization::make_nvp("control_points", control_points_);
  }
};  // End struct bezier_curve

template <typename T, typename N, bool S, typename P>
bezier_curve<T, N, S, P> operator+(const bezier_curve<T, N, S, P>& p1,
                                   const bezier_curve<T, N, S, P>& p2) {
  bezier_curve<T, N, S, P> res(p1);
  return res += p2;
}

template <typename T, typename N, bool S, typename P>
bezier_curve<T, N, S, P> operator-(const bezier_curve<T, N, S, P>& p1) {
  std::vector<typename bezier_curve<T, N, S, P>::point_t> ts;
  for (std::size_t i = 0; i <= p1.degree(); ++i) {
    ts.push_back(bezier_curve<T, N, S, P>::point_t::Zero(p1.dim()));
  }
  bezier_curve<T, N, S, P> res(ts.begin(), ts.end(), p1.min(), p1.max());
  res -= p1;
  return res;
}

template <typename T, typename N, bool S, typename P>
bezier_curve<T, N, S, P> operator-(const bezier_curve<T, N, S, P>& p1,
                                   const bezier_curve<T, N, S, P>& p2) {
  bezier_curve<T, N, S, P> res(p1);
  return res -= p2;
}

template <typename T, typename N, bool S, typename P>
bezier_curve<T, N, S, P> operator-(
    const bezier_curve<T, N, S, P>& p1,
    const typename bezier_curve<T, N, S, P>::point_t& point) {
  bezier_curve<T, N, S, P> res(p1);
  return res -= point;
}

template <typename T, typename N, bool S, typename P>
bezier_curve<T, N, S, P> operator-(
    const typename bezier_curve<T, N, S, P>::point_t& point,
    const bezier_curve<T, N, S, P>& p1) {
  bezier_curve<T, N, S, P> res(-p1);
  return res += point;
}

template <typename T, typename N, bool S, typename P>
bezier_curve<T, N, S, P> operator+(
    const bezier_curve<T, N, S, P>& p1,
    const typename bezier_curve<T, N, S, P>::point_t& point) {
  bezier_curve<T, N, S, P> res(p1);
  return res += point;
}

template <typename T, typename N, bool S, typename P>
bezier_curve<T, N, S, P> operator+(
    const typename bezier_curve<T, N, S, P>::point_t& point,
    const bezier_curve<T, N, S, P>& p1) {
  bezier_curve<T, N, S, P> res(p1);
  return res += point;
}

template <typename T, typename N, bool S, typename P>
bezier_curve<T, N, S, P> operator/(const bezier_curve<T, N, S, P>& p1,
                                   const double k) {
  bezier_curve<T, N, S, P> res(p1);
  return res /= k;
}

template <typename T, typename N, bool S, typename P>
bezier_curve<T, N, S, P> operator*(const bezier_curve<T, N, S, P>& p1,
                                   const double k) {
  bezier_curve<T, N, S, P> res(p1);
  return res *= k;
}

template <typename T, typename N, bool S, typename P>
bezier_curve<T, N, S, P> operator*(const double k,
                                   const bezier_curve<T, N, S, P>& p1) {
  bezier_curve<T, N, S, P> res(p1);
  return res *= k;
}

}  // namespace ndcurves

DEFINE_CLASS_TEMPLATE_VERSION(
    SINGLE_ARG(typename Time, typename Numeric, bool Safe, typename Point),
    SINGLE_ARG(ndcurves::bezier_curve<Time, Numeric, Safe, Point>))

#endif  //_CLASS_BEZIERCURVE
