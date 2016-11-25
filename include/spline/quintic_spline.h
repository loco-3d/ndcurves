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

#include "spline_curve.h"

#include <stdexcept>

namespace spline
{
/// \class quintic_spline
/// \brief Represents a quintic spline defined on the interval
/// [tBegin, tEnd]. It follows the equation
/// x(t) = a + b(t - t_min_) + c(t - t_min_)^2 + d(t - t_min_)^3 + e(t - t_min_)^4  + f(t - t_min_)^5
///
template<typename Time= double, typename Numeric=Time, std::size_t Dim=3, bool Safe=false
, typename Point= Eigen::Matrix<Numeric, Dim, 1>, typename T_Point =std::vector<Point,Eigen::aligned_allocator<Point> > >
struct quintic_spline : public spline_curve<Time, Numeric, Dim,5, Safe, Point, T_Point>
{
    typedef Point 	point_t;
    typedef T_Point t_point_t;
    typedef Time 	time_t;
    typedef Numeric	num_t;
    typedef spline_curve<Time, Numeric, Dim,5, Safe, Point, T_Point> spline_curve_t;
/* Constructors - destructors */
    public:
    ///\brief Constructor
    quintic_spline(point_t const& a, point_t const& b, point_t const& c, point_t const &d, point_t const &e,  point_t const &f, time_t min, time_t max)
        :spline_curve_t(makeVector(a,b,c,d,e,f), min, max) {}


    ///\brief Constructor
    quintic_spline(const T_Point& coefficients, time_t min, time_t max)
        :spline_curve_t(coefficients, min, max) {}


    ///\brief Constructor
    template<typename In>
    quintic_spline(In zeroOrderCoefficient, In out, time_t min, time_t max)
        :spline_curve_t(zeroOrderCoefficient, out, min, max) {}

    ///\brief Destructor
    ~quintic_spline()
    {
        // NOTHING
    }

    private:
    //quintic_spline(const quintic_spline&);
    quintic_spline& operator=(const quintic_spline&);
/* Constructors - destructors */

    /*Operations*/
    T_Point makeVector(point_t const& a, point_t const& b, point_t const& c,
                       point_t const &d, point_t const& e, point_t const& f)
    {
        T_Point res;
        res.push_back(a);res.push_back(b);res.push_back(c);
        res.push_back(d);res.push_back(e);res.push_back(f);
        return res;
    }

    /*Operations*/
    }; //class quintic_spline
}
#endif //_STRUCT_QUINTIC_SPLINE

