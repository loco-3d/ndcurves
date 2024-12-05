/**
 * \file exact_cubic.h
 * \brief class allowing to create an Exact cubic spline.
 * \author Steve T.
 * \version 0.1
 * \date 06/17/2013
 *
 * This file contains definitions for the ExactCubic class.
 * Given a set of waypoints (x_i*) and timestep (t_i), it provides the unique
 * set of cubic splines fulfulling those 4 restrictions :
 * - x_i(t_i) = x_i* ; this means that the curve passes trough each waypoint
 * - x_i(t_i+1) = x_i+1* ;
 * - its derivative is continous at t_i+1
 * - its acceleration is continous at t_i+1
 * more details in paper "Task-Space Trajectories via Cubic Spline Optimization"
 * By J. Zico Kolter and Andrew Y.ng (ICRA 2009)
 */

#ifndef _CLASS_EXACTCUBIC
#define _CLASS_EXACTCUBIC

#include <functional>
#include <vector>

#include "MathDefs.h"
#include "curve_abc.h"
#include "curve_constraint.h"
#include "piecewise_curve.h"
#include "polynomial.h"

namespace ndcurves {
/// \class ExactCubic.
/// \brief Represents a set of cubic splines defining a continuous function
/// crossing each of the waypoint given in its initialization.
///
template <
    typename Time = double, typename Numeric = Time, bool Safe = false,
    typename Point = Eigen::Matrix<Numeric, Eigen::Dynamic, 1>,
    typename T_Point = std::vector<Point, Eigen::aligned_allocator<Point> >,
    typename SplineBase = polynomial<Time, Numeric, Safe, Point, T_Point> >
struct exact_cubic : public piecewise_curve<Time, Numeric, Safe, Point> {
  typedef Point point_t;
  typedef const Eigen::Ref<const point_t> point_ref_t;
  typedef T_Point t_point_t;
  typedef Eigen::Matrix<Numeric, Eigen::Dynamic, Eigen::Dynamic> MatrixX;
  typedef Eigen::Matrix<Numeric, 3, 3> Matrix3;
  typedef Time time_t;
  typedef Numeric num_t;
  typedef SplineBase spline_t;
  typedef typename std::vector<spline_t> t_spline_t;
  typedef typename t_spline_t::iterator it_spline_t;
  typedef typename t_spline_t::const_iterator cit_spline_t;
  typedef curve_constraints<Point> spline_constraints;

  typedef exact_cubic<Time, Numeric, Safe, point_t, T_Point, SplineBase>
      exact_cubic_t;
  typedef curve_abc<Time, Numeric, Safe, point_t> curve_abc_t;  // parent class
  typedef piecewise_curve<Time, Numeric, Safe, point_t> piecewise_curve_t;
  typedef polynomial<Time, Numeric, Safe, point_t> polynomial_t;
  typedef typename piecewise_curve_t::t_curve_ptr_t t_curve_ptr_t;

  /* Constructors - destructors */
 public:
  /// \brief Empty constructor. Add at least one curve to call other class
  /// functions.
  ///
  exact_cubic() : piecewise_curve_t() {}

  /// \brief Constructor.
  /// \param wayPointsBegin : an iterator pointing to the first element of a
  /// waypoint container. \param wayPointsEns   : an iterator pointing to the
  /// last element of a waypoint container.
  ///
  template <typename In>
  exact_cubic(In wayPointsBegin, In wayPointsEnd) : piecewise_curve_t() {
    t_spline_t subSplines = computeWayPoints<In>(wayPointsBegin, wayPointsEnd);
    for (cit_spline_t it = subSplines.begin(); it != subSplines.end(); ++it) {
      this->add_curve(*it);
    }
  }

  /// \brief Constructor.
  /// \param wayPointsBegin : an iterator pointing to the first element of a
  /// waypoint container. \param wayPointsEns   : an iterator pointing to the
  /// last element of a waypoint container. \param constraints    : constraints
  /// on the init and end velocity / accelerations of the spline.
  ///
  template <typename In>
  exact_cubic(In wayPointsBegin, In wayPointsEnd,
              const spline_constraints& constraints)
      : piecewise_curve_t() {
    t_spline_t subSplines =
        computeWayPoints<In>(wayPointsBegin, wayPointsEnd, constraints);
    for (cit_spline_t it = subSplines.begin(); it != subSplines.end(); ++it) {
      this->add_curve(*it);
    }
  }

  /// \brief Constructor.
  /// \param subSplines: vector of subSplines.
  exact_cubic(const t_spline_t& subSplines) : piecewise_curve_t() {
    for (cit_spline_t it = subSplines.begin(); it != subSplines.end(); ++it) {
      this->add_curve(*it);
    }
  }

  exact_cubic(const t_curve_ptr_t& subSplines)
      : piecewise_curve_t(subSplines) {}

  /// \brief Copy Constructor.
  exact_cubic(const exact_cubic& other) : piecewise_curve_t(other) {}

  /// \brief Destructor.
  virtual ~exact_cubic() {}

  std::size_t getNumberSplines() { return this->getNumberCurves(); }

  spline_t getSplineAt(std::size_t index) {
    std::shared_ptr<spline_t> s_ptr =
        std::dynamic_pointer_cast<spline_t>(this->curves_.at(index));
    if (s_ptr)
      return *s_ptr;
    else
      throw std::runtime_error(
          "Parent piecewise curve do not contain only curves created from "
          "exact_cubic class methods");
  }

 private:
  static polynomial_t create_cubic(point_ref_t a, point_ref_t b, point_ref_t c,
                                   point_ref_t d, const time_t t_min,
                                   const time_t t_max) {
    typename polynomial_t::t_point_t coeffs;
    coeffs.push_back(a);
    coeffs.push_back(b);
    coeffs.push_back(c);
    coeffs.push_back(d);
    return polynomial_t(coeffs.begin(), coeffs.end(), t_min, t_max);
  }
  static polynomial_t create_quintic(point_ref_t a, point_ref_t b,
                                     point_ref_t c, point_ref_t d,
                                     point_ref_t e, point_ref_t f,
                                     const time_t t_min, const time_t t_max) {
    typename polynomial_t::t_point_t coeffs;
    coeffs.push_back(a);
    coeffs.push_back(b);
    coeffs.push_back(c);
    coeffs.push_back(d);
    coeffs.push_back(e);
    coeffs.push_back(f);
    return polynomial_t(coeffs.begin(), coeffs.end(), t_min, t_max);
  }

  /// \brief Compute polynom of exact cubic spline from waypoints.
  /// Compute the coefficients of polynom as in paper : "Task-Space Trajectories
  /// via Cubic Spline Optimization".<br>
  /// \f$x_i(t)=a_i+b_i(t-t_i)+c_i(t-t_i)^2\f$<br>
  /// with \f$a=x\f$, \f$H_1b=H_2x\f$, \f$c=H_3x+H_4b\f$, \f$d=H_5x+H_6b\f$.<br>
  /// The matrices \f$H\f$ are defined as in the paper in Appendix A.
  ///
  template <typename In>
  t_spline_t computeWayPoints(In wayPointsBegin, In wayPointsEnd) const {
    const std::size_t dim = wayPointsBegin->second.size();
    const std::size_t size = std::distance(wayPointsBegin, wayPointsEnd);
    if (Safe && size < 1) {
      throw std::length_error(
          "size of waypoints must be superior to 0");  // TODO
    }
    t_spline_t subSplines;
    subSplines.reserve(size);
    // refer to the paper to understand all this.
    MatrixX h1 = MatrixX::Zero(size, size);
    MatrixX h2 = MatrixX::Zero(size, size);
    MatrixX h3 = MatrixX::Zero(size, size);
    MatrixX h4 = MatrixX::Zero(size, size);
    MatrixX h5 = MatrixX::Zero(size, size);
    MatrixX h6 = MatrixX::Zero(size, size);
    MatrixX a = MatrixX::Zero(size, dim);
    MatrixX b = MatrixX::Zero(size, dim);
    MatrixX c = MatrixX::Zero(size, dim);
    MatrixX d = MatrixX::Zero(size, dim);
    MatrixX x = MatrixX::Zero(size, dim);
    In it(wayPointsBegin), next(wayPointsBegin);
    ++next;
    // Fill the matrices H as specified in the paper.
    for (std::size_t i(0); next != wayPointsEnd; ++next, ++it, ++i) {
      num_t const dTi((*next).first - (*it).first);
      num_t const dTi_sqr(dTi * dTi);
      num_t const dTi_cube(dTi_sqr * dTi);
      // filling matrices values
      h3(i, i) = -3 / dTi_sqr;
      h3(i, i + 1) = 3 / dTi_sqr;
      h4(i, i) = -2 / dTi;
      h4(i, i + 1) = -1 / dTi;
      h5(i, i) = 2 / dTi_cube;
      h5(i, i + 1) = -2 / dTi_cube;
      h6(i, i) = 1 / dTi_sqr;
      h6(i, i + 1) = 1 / dTi_sqr;
      if (i + 2 < size) {
        In it2(next);
        ++it2;
        num_t const dTi_1((*it2).first - (*next).first);
        num_t const dTi_1sqr(dTi_1 * dTi_1);
        // this can be optimized but let's focus on clarity as long as not
        // needed
        h1(i + 1, i) = 2 / dTi;
        h1(i + 1, i + 1) = 4 / dTi + 4 / dTi_1;
        h1(i + 1, i + 2) = 2 / dTi_1;
        h2(i + 1, i) = -6 / dTi_sqr;
        h2(i + 1, i + 1) = (6 / dTi_1sqr) - (6 / dTi_sqr);
        h2(i + 1, i + 2) = 6 / dTi_1sqr;
      }
      x.row(i) = (*it).second.transpose();
    }
    // adding last x
    x.row(size - 1) = (*it).second.transpose();
    // Compute coefficients of polynom.
    a = x;
    PseudoInverse(h1);
    b = h1 * h2 * x;  // h1 * b = h2 * x => b = (h1)^-1 * h2 * x
    c = h3 * x + h4 * b;
    d = h5 * x + h6 * b;
    // create splines along waypoints.
    it = wayPointsBegin, next = wayPointsBegin;
    ++next;
    for (int i = 0; next != wayPointsEnd; ++i, ++it, ++next) {
      subSplines.push_back(create_cubic(a.row(i), b.row(i), c.row(i), d.row(i),
                                        (*it).first, (*next).first));
    }
    return subSplines;
  }

  template <typename In>
  t_spline_t computeWayPoints(In wayPointsBegin, In wayPointsEnd,
                              const spline_constraints& constraints) const {
    std::size_t const size(std::distance(wayPointsBegin, wayPointsEnd));
    if (Safe && size < 1) {
      throw std::length_error(
          "number of waypoints should be superior to one");  // TODO
    }
    t_spline_t subSplines;
    subSplines.reserve(size - 1);
    spline_constraints cons = constraints;
    In it(wayPointsBegin), next(wayPointsBegin), end(wayPointsEnd - 1);
    ++next;
    for (std::size_t i(0); next != end; ++next, ++it, ++i) {
      compute_one_spline<In>(it, next, cons, subSplines);
    }
    compute_end_spline<In>(it, next, cons, subSplines);
    return subSplines;
  }

  template <typename In>
  void compute_one_spline(In wayPointsBegin, In wayPointsNext,
                          spline_constraints& constraints,
                          t_spline_t& subSplines) const {
    const point_t &a0 = wayPointsBegin->second, a1 = wayPointsNext->second;
    const point_t &b0 = constraints.init_vel, c0 = constraints.init_acc / 2.;
    const num_t &init_t = wayPointsBegin->first, end_t = wayPointsNext->first;
    const num_t dt = end_t - init_t, dt_2 = dt * dt, dt_3 = dt_2 * dt;
    const point_t d0 = (a1 - a0 - b0 * dt - c0 * dt_2) / dt_3;
    subSplines.push_back(create_cubic(a0, b0, c0, d0, init_t, end_t));
    constraints.init_vel = subSplines.back().derivate(end_t, 1);
    constraints.init_acc = subSplines.back().derivate(end_t, 2);
  }

  template <typename In>
  void compute_end_spline(In wayPointsBegin, In wayPointsNext,
                          spline_constraints& constraints,
                          t_spline_t& subSplines) const {
    const std::size_t dim = wayPointsBegin->second.size();
    const point_t &a0 = wayPointsBegin->second, a1 = wayPointsNext->second;
    const point_t &b0 = constraints.init_vel, b1 = constraints.end_vel,
                  c0 = constraints.init_acc / 2., c1 = constraints.end_acc;
    const num_t &init_t = wayPointsBegin->first, end_t = wayPointsNext->first;
    const num_t dt = end_t - init_t, dt_2 = dt * dt, dt_3 = dt_2 * dt,
                dt_4 = dt_3 * dt, dt_5 = dt_4 * dt;
    // solving a system of four linear eq with 4 unknows: d0, e0
    const point_t alpha_0 = a1 - a0 - b0 * dt - c0 * dt_2;
    const point_t alpha_1 = b1 - b0 - 2 * c0 * dt;
    const point_t alpha_2 = c1 - 2 * c0;
    const num_t x_d_0 = dt_3, x_d_1 = 3 * dt_2, x_d_2 = 6 * dt;
    const num_t x_e_0 = dt_4, x_e_1 = 4 * dt_3, x_e_2 = 12 * dt_2;
    const num_t x_f_0 = dt_5, x_f_1 = 5 * dt_4, x_f_2 = 20 * dt_3;
    point_t d, e, f;
    MatrixX rhs = MatrixX::Zero(3, dim);
    rhs.row(0) = alpha_0;
    rhs.row(1) = alpha_1;
    rhs.row(2) = alpha_2;
    Matrix3 eq = Matrix3::Zero(3, 3);
    eq(0, 0) = x_d_0;
    eq(0, 1) = x_e_0;
    eq(0, 2) = x_f_0;
    eq(1, 0) = x_d_1;
    eq(1, 1) = x_e_1;
    eq(1, 2) = x_f_1;
    eq(2, 0) = x_d_2;
    eq(2, 1) = x_e_2;
    eq(2, 2) = x_f_2;
    rhs = eq.inverse().eval() * rhs;
    d = rhs.row(0);
    e = rhs.row(1);
    f = rhs.row(2);
    subSplines.push_back(create_quintic(a0, b0, c0, d, e, f, init_t, end_t));
  }

 public:
  // Serialization of the class
  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version) {
    if (version) {
      // Do something depending on version ?
    }
    ar& BOOST_SERIALIZATION_BASE_OBJECT_NVP(piecewise_curve_t);
  }
};
}  // namespace ndcurves

DEFINE_CLASS_TEMPLATE_VERSION(
    SINGLE_ARG(typename Time, typename Numeric, bool Safe, typename Point,
               typename T_Point, typename SplineBase),
    SINGLE_ARG(
        ndcurves::exact_cubic<Time, Numeric, Safe, Point, T_Point, SplineBase>))
#endif  //_CLASS_EXACTCUBIC
