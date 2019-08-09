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

#include "polynomial.h"

#include <stdexcept>

namespace curves
{
	/// \brief Creates coefficient vector of a cubic spline defined on the interval
	/// \f$[t_{min}, t_{max}]\f$. It follows the equation : <br>
	/// \f$ x(t) = a + b(t - t_{min}) + c(t - t_{min})^2 + d(t - t_{min})^3 \f$ where \f$ t \in [t_{min}, t_{max}] \f$
	/// with a, b, c and d the control points.
	///
	template<typename Point, typename T_Point>
	T_Point make_cubic_vector(Point const& a, Point const& b, Point const& c, Point const &d)
	{
		T_Point res;
		res.push_back(a);res.push_back(b);res.push_back(c);res.push_back(d);
		return res;
	}

	template<typename Time, typename Numeric, std::size_t Dim, bool Safe, 
					 typename Point, typename T_Point>
	polynomial<Time,Numeric,Dim,Safe,Point,T_Point> create_cubic(Point const& a, Point const& b, Point const& c, Point const &d,
																															 const Time t_min, const Time t_max)
	{
		T_Point coeffs = make_cubic_vector<Point, T_Point>(a,b,c,d);
		return polynomial<Time,Numeric,Dim,Safe,Point,T_Point>(coeffs.begin(),coeffs.end(), t_min, t_max);
	}
} // namespace curves
#endif //_STRUCT_CUBICSPLINE

