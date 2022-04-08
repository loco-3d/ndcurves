#ifndef _CLASS_CURVE_CONVERSION
#define _CLASS_CURVE_CONVERSION

#include <iostream>
#include <stdexcept>
#include <vector>

#include "MathDefs.h"
#include "bernstein.h"
#include "curve_abc.h"
#include "curve_constraint.h"

namespace ndcurves {
/// \brief Converts a cubic hermite spline or a bezier curve to a polynomial.
/// \param curve   : the bezier curve/cubic hermite spline defined between
/// [Tmin,Tmax] to convert. \return the equivalent polynomial.
template <typename Polynomial>
Polynomial polynomial_from_curve(
    const typename Polynomial::curve_abc_t& curve) {
  typedef typename Polynomial::t_point_t t_point_t;
  typedef typename Polynomial::num_t num_t;
  t_point_t coefficients;
  coefficients.push_back(curve(curve.min()));
  num_t fact = 1;
  for (std::size_t i = 1; i <= curve.degree(); ++i) {
    fact *= (num_t)i;
    coefficients.push_back(curve.derivate(curve.min(), i) / fact);
  }

  return Polynomial(coefficients, curve.min(), curve.max());
}

/// \brief Converts a cubic hermite spline or polynomial of order 3 or less to a
/// cubic bezier curve. \param curve   : the polynomial of order 3 or less/cubic
/// hermite spline defined between [Tmin,Tmax] to convert. \return the
/// equivalent cubic bezier curve.
template <typename Bezier>
Bezier bezier_from_curve(const typename Bezier::curve_abc_t& curve) {
  if (curve.degree() > 3)
    throw std::invalid_argument(
        "bezier_from_curve is only implemented for curves of degree <= 3.");
  typedef typename Bezier::point_t point_t;
  typedef typename Bezier::t_point_t t_point_t;
  typedef typename Bezier::num_t num_t;
  num_t T_min = curve.min();
  num_t T_max = curve.max();
  num_t T = T_max - T_min;
  // Positions and derivatives
  point_t p0 = curve(T_min);
  point_t p1 = curve(T_max);
  point_t m0 = curve.derivate(T_min, 1);
  point_t m1 = curve.derivate(T_max, 1);
  // Convert to bezier control points
  // for t in [Tmin,Tmax] and T=Tmax-Tmin : x'(0)=3(b_p1-b_p0)/T and
  // x'(1)=3(b_p3-b_p2)/T so : m0=3(b_p1-b_p0)/T and m1=3(b_p3-b_p2)/T
  // <=> b_p1=T(m0/3)+b_p0 and b_p2=-T(m1/3)+b_p3
  point_t b_p0 = p0;
  point_t b_p3 = p1;
  point_t b_p1 = T * m0 / 3 + b_p0;
  point_t b_p2 = -T * m1 / 3 + b_p3;
  t_point_t control_points;
  control_points.push_back(b_p0);
  control_points.push_back(b_p1);
  control_points.push_back(b_p2);
  control_points.push_back(b_p3);
  return Bezier(control_points.begin(), control_points.end(), curve.min(),
                curve.max());
}

/// \brief Converts a polynomial of order 3 or less/cubic bezier curve to a
/// cubic hermite spline. \param curve   : the polynomial of order 3 or
/// less/cubic bezier curve defined between [Tmin,Tmax] to convert. \return the
/// equivalent cubic hermite spline.
template <typename Hermite>
Hermite hermite_from_curve(const typename Hermite::curve_abc_t& curve) {
  if (curve.degree() > 3)
    throw std::invalid_argument(
        "hermite_from_curve is only implemented for curves of degree <= 3.");
  typedef typename Hermite::pair_point_tangent_t pair_point_tangent_t;
  typedef typename Hermite::t_pair_point_tangent_t t_pair_point_tangent_t;
  typedef typename Hermite::point_t point_t;
  typedef typename Hermite::num_t num_t;
  num_t T_min = curve.min();
  num_t T_max = curve.max();
  // Positions and derivatives
  point_t p0 = curve(T_min);
  point_t p1 = curve(T_max);
  point_t m0 = curve.derivate(T_min, 1);
  point_t m1 = curve.derivate(T_max, 1);
  // Create pairs pos/vel
  pair_point_tangent_t pair0(p0, m0);
  pair_point_tangent_t pair1(p1, m1);
  t_pair_point_tangent_t control_points;
  control_points.push_back(pair0);
  control_points.push_back(pair1);
  std::vector<double> time_control_points;
  time_control_points.push_back(T_min);
  time_control_points.push_back(T_max);
  return Hermite(control_points.begin(), control_points.end(),
                 time_control_points);
}
}  // namespace ndcurves
#endif  //_CLASS_CURVE_CONVERSION
