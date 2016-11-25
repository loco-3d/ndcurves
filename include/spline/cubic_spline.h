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


#ifndef _STRUCT_CUBICSPLINE
#define _STRUCT_CUBICSPLINE

#include "MathDefs.h"

#include "spline_curve.h"

#include <stdexcept>

namespace spline
{
/// \brief Creates coefficient vector of a cubic spline defined on the interval
/// [tBegin, tEnd]. It follows the equation
/// x(t) = a + b(t - t_min_) + c(t - t_min_)^2 + d(t - t_min_)^3
///
template<typename Point, typename T_Point>
T_Point make_cubic_vector(Point const& a, Point const& b, Point const& c, Point const &d)
{
    T_Point res;
    res.push_back(a);res.push_back(b);res.push_back(c);res.push_back(d);
    return res;
}
}
#endif //_STRUCT_CUBICSPLINE

