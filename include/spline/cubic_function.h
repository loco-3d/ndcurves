/**
* \file cubic_function.h
* \brief Definition of a cubic spline.
* \author Steve T.
* \version 0.1
* \date 06/17/2013
*
* This file contains definitions for the CubicFunction struct.
* It allows the creation and evaluation of natural
* smooth cubic splines of arbitrary dimension
*/


#ifndef _STRUCT_CUBICFUNCTION
#define _STRUCT_CUBICFUNCTION

#include "MathDefs.h"

#include "spline_curve.h"

#include <stdexcept>

namespace spline
{
/// \class CubicFunction
/// \brief Represents a cubic spline defined on the interval 
/// [tBegin, tEnd]. It follows the equation
/// x(t) = a + b(t - t_min_) + c(t - t_min_)^2 + d(t - t_min_)^3 
///
template<typename Time= double, typename Numeric=Time, std::size_t Dim=3, bool Safe=false
, typename Point= Eigen::Matrix<Numeric, Dim, 1>, typename T_Point =std::vector<Point,Eigen::aligned_allocator<Point> > >
struct cubic_function : public spline_curve<Time, Numeric, Dim,3, Safe, Point, T_Point>
{
    typedef Point 	point_t;
    typedef T_Point t_point_t;
  	typedef Time 	time_t;
    typedef Numeric	num_t;
    typedef spline_curve<Time, Numeric, Dim,3, Safe, Point, T_Point> spline_curve_t;
/* Constructors - destructors */
	public:
	///\brief Constructor
	cubic_function(point_t const& a, point_t const& b, point_t const& c, point_t const &d, time_t min, time_t max)
        :spline_curve_t(makeVector(a,b,c,d), min, max) {}


    ///\brief Constructor
    cubic_function(const T_Point& coefficients, time_t min, time_t max)
        :spline_curve_t(coefficients, min, max) {}


    ///\brief Constructor
    template<typename In>
    cubic_function(In zeroOrderCoefficient, In out, time_t min, time_t max)
        :spline_curve_t(zeroOrderCoefficient, out, min, max) {}

    ///\brief Destructor
    ~cubic_function()
    {
        // NOTHING
    }

    private:
    cubic_function(const cubic_function&);
    cubic_function& operator=(const cubic_function&);
/* Constructors - destructors */

    /*Operations*/
    T_Point makeVector(point_t const& a, point_t const& b, point_t const& c, point_t const &d)
    {
        T_Point res;
        res.push_back(a);res.push_back(b);res.push_back(c);res.push_back(d);
        return res;
    }

    /*Operations*/
    }; //class CubicFunction
}
#endif //_STRUCT_CUBICFUNCTION

