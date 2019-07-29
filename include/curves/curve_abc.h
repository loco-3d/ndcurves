/**
* \file curve_abc.h
* \brief interface for a Curve of arbitrary dimension.
* \author Steve T.
* \version 0.1
* \date 06/17/2013
*
* Interface for a curve
*/


#ifndef _STRUCT_CURVE_ABC
#define _STRUCT_CURVE_ABC

#include "MathDefs.h"

#include <functional>

namespace curves
{
/// \struct curve_abc.
/// \brief Represents a curve of dimension Dim.
/// If value of parameter Safe is false, no verification is made on the evaluation of the curve.
template<typename Time= double, typename Numeric=Time, bool Safe=false
, typename Point= Eigen::Matrix<Numeric, Eigen::Dynamic, 1> >
struct  curve_abc : std::unary_function<Time, Point>
{
    typedef Point   point_t;
    typedef Time    time_t;

/* Constructors - destructors */
    public:
    /// \brief Constructor.
    curve_abc(){}

    /// \brief Destructor.
    virtual ~curve_abc(){}
/* Constructors - destructors */

/*Operations*/
    public:
    ///  \brief Evaluation of the cubic spline at time t.
    ///  \param t : time when to evaluate the spine
    ///  \return \f$x(t)\f$, point corresponding on curve at time t.
    virtual point_t operator()(const time_t t) const = 0;


    /// \brief Evaluate the derivative of order N of curve at time t.
    /// \param t : time when to evaluate the spline.
    /// \param order : order of derivative.
    /// \return \f$\frac{d^Nx(t)}{dt^N}\f$, point corresponding on derivative curve of order N at time t.
    virtual point_t derivate(const time_t t, const std::size_t order) const = 0;
/*Operations*/

/*Helpers*/
    public:
    /// \brief Get the minimum time for which the curve is defined.
    /// \return \f$t_{min}\f$, lower bound of time range.
    virtual time_t min() const = 0;
    /// \brief Get the maximum time for which the curve is defined.
    /// \return \f$t_{max}\f$, upper bound of time range.
    virtual time_t max() const = 0;

    std::pair<time_t, time_t> timeRange() {return std::make_pair(min(), max());}
/*Helpers*/

    };
} // namespace curves
#endif //_STRUCT_CURVE_ABC

