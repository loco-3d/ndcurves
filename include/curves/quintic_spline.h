/**
 * \file cubic_spline.h
 * \brief Definition of a cubic spline.
 * \author Steve T.
 * \version 0.1
 * \date 06/17/2013
 *
 * This file contains definitions for the CubicFunction struct.
 * It allows the creation and evaluation of natural
 * smooth cubic splines of arbitrary dimension
 */

#ifndef _STRUCT_QUINTIC_SPLINE
#define _STRUCT_QUINTIC_SPLINE

#include "MathDefs.h"

#include "polynomial.h"

#include <stdexcept>

namespace curves {
/// \brief Creates coefficient vector of a quintic spline defined on the interval
/// \f$[t_{min}, t_{max}]\f$. It follows the equation :<br>
/// \f$ x(t) = a + b(t - t_{min}) + c(t - t_{min})^2 + d(t - t_{min})^3 + e(t - t_{min})^4  + f(t - t_{min})^5 \f$ <br>
/// where \f$ t \in [t_{min}, t_{max}] \f$.
///
template <typename Point, typename T_Point>
T_Point make_quintic_vector(Point const& a, Point const& b, Point const& c, Point const& d, Point const& e,
                            Point const& f) {
  T_Point res;
  res.push_back(a);
  res.push_back(b);
  res.push_back(c);
  res.push_back(d);
  res.push_back(e);
  res.push_back(f);
  return res;
}

template <typename Time, typename Numeric, bool Safe, typename Point, typename T_Point>
polynomial<Time, Numeric, Safe, Point, T_Point> create_quintic(Point const& a, Point const& b, Point const& c,
                                                               Point const& d, Point const& e, Point const& f,
                                                               const Time t_min, const Time t_max) {
  T_Point coeffs = make_quintic_vector<Point, T_Point>(a, b, c, d, e, f);
  return polynomial<Time, Numeric, Safe, Point, T_Point>(coeffs.begin(), coeffs.end(), t_min, t_max);
}
}  // namespace curves
#endif  //_STRUCT_QUINTIC_SPLINE
