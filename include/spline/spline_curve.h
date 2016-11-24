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


#ifndef _STRUCT_SPLINE
#define _STRUCT_SPLINE

#include "MathDefs.h"

#include "curve_abc.h"

#include <iostream>
#include <algorithm>
#include <functional>
#include <stdexcept>

namespace spline
{
/// \class spline_curve
/// \brief Represents a spline curve of arbitrary order defined on the interval
/// [tBegin, tEnd]. It follows the equation
/// x(t) = a + b(t - t_min_) + ... + d(t - t_min_)^N, where N is the order
///
template<typename Time= double, typename Numeric=Time, std::size_t Dim=3, std::size_t Order=3, bool Safe=false,
         typename Point= Eigen::Matrix<Numeric, Dim, 1>, typename T_Point =std::vector<Point,Eigen::aligned_allocator<Point> > >
struct spline_curve : public curve_abc<Time, Numeric, Dim, Safe, Point>
{
    typedef Point 	point_t;
    typedef T_Point  t_point_t;
    typedef typename t_point_t::const_iterator cit_point_t;
    typedef Time 	time_t;
    typedef Numeric	num_t;
/* Constructors - destructors */
    public:
    ///\brief Constructor
    ///\param coefficients : a container containing all coefficients of the spline, starting
    /// with the zero order coefficient, up to the highest order
    ///\param min: LOWER bound on interval definition of the spline
    ///\param max: UPPER bound on interval definition of the spline
    spline_curve(const T_Point& coefficients, time_t min, time_t max)
        :coefficients_(coefficients), t_min_(min), t_max_(max), dim_(Dim), order_(Order)
    {
        if(t_min_ > t_max_ && Safe)
        {
            std::out_of_range("TODO");
        }
    }

    ///\brief Constructor
    ///\param zeroOrderCoefficient : an iterator pointing to the first element of a structure containing the coefficients
    /// it corresponds to the zero degree coefficient
    ///\param out   : an iterator pointing to the last element of a structure ofcoefficients
    ///\param min: LOWER bound on interval definition of the spline
    ///\param max: UPPER bound on interval definition of the spline
    template<typename In>
    spline_curve(In zeroOrderCoefficient, In out, time_t min, time_t max)
        :coefficients_(init_coeffs(zeroOrderCoefficient, out)), t_min_(min), t_max_(max), dim_(Dim), order_(Order)
    {
        if(t_min_ > t_max_ && Safe)
        {
            std::out_of_range("TODO");
        }
    }

    ///\brief Destructor
    ~spline_curve()
    {
        // NOTHING
    }

    private:
    spline_curve(const spline_curve&);
    spline_curve& operator=(const spline_curve&);
/* Constructors - destructors */

/*Operations*/
    public:
    ///  \brief Evaluation of the cubic spline at time t.
    ///  \param t : the time when to evaluate the spine
    ///  \param return : the value x(t)
    virtual point_t operator()(time_t t) const
    {
        if((t < t_min_ || t > t_max_) && Safe){ throw std::out_of_range("TODO");}
        time_t const dt (t-t_min_);
        time_t cdt(dt);
        cit_point_t cit = coefficients_.begin();
        point_t currentPoint_ = *cit; ++cit;
        for(; cit != coefficients_.end(); ++cit, cdt*=dt)
            currentPoint_ += cdt *(*cit);
        return currentPoint_;
    }
/*Operations*/

/*Helpers*/
    public:
    ///  \brief Returns the minimum time for wich curve is defined
    num_t virtual min() const {return t_min_;}
    ///  \brief Returns the maximum time for wich curve is defined
    num_t virtual max() const {return t_max_;}
/*Helpers*/

/*Attributes*/
    public:
    const t_point_t coefficients_;
    const time_t t_min_, t_max_;
    const std::size_t dim_;
    const std::size_t order_;
/*Attributes*/

    private:
    template<typename In>
    t_point_t init_coeffs(In zeroOrderCoefficient, In highestOrderCoefficient)
    {
        t_point_t res(Order+1);
        std::copy(zeroOrderCoefficient, highestOrderCoefficient, res.begin());
        return res;
    }
}; //class spline_curve
}
#endif //_STRUCT_SPLINE

