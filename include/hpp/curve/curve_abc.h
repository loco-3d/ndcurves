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

namespace curve
{
/// \struct curve_abc
/// \brief Represents a curve of dimension Dim
/// is Safe is false, no verification is made on the evaluation of the curve.
template<typename Time= double, typename Numeric=Time, std::size_t Dim=3, bool Safe=false
, typename Point= Eigen::Matrix<Numeric, Dim, 1> >
struct  curve_abc : std::unary_function<Time, Point>
{
	typedef Point 	point_t;
	typedef Time 	time_t;

/* Constructors - destructors */
	public:
	///\brief Constructor
    curve_abc(){}

	///\brief Destructor
    virtual ~curve_abc(){}
/* Constructors - destructors */

/*Operations*/
	public:
	///  \brief Evaluation of the cubic spline at time t.
	///  \param t : the time when to evaluate the spine
	///  \param return : the value x(t)
    virtual point_t operator()(const time_t t) const = 0;


    ///  \brief Evaluation of the derivative spline at time t.
    ///  \param t : the time when to evaluate the spline
    ///  \param order : order of the derivative
    ///  \param return : the value x(t)
    virtual point_t derivate(const time_t t, const std::size_t order) const = 0;
/*Operations*/

/*Helpers*/
	public:
	///  \brief Returns the minimum time for wich curve is defined
	virtual time_t min() const = 0;
	///  \brief Returns the maximum time for wich curve is defined
	virtual time_t max() const = 0;

    std::pair<time_t, time_t> timeRange() {return std::make_pair(min(), max());}
/*Helpers*/

	};
}
#endif //_STRUCT_CURVE_ABC

