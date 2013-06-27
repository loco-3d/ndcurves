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

#include "Exports.h"
#include "MathDefs.h"

#include <functional>

namespace spline
{
/// \struct curve_abc
/// \brief Represents a curve of dimension Dim
/// is Safe is false, no verification is made on the evaluation of the curve.
template<typename Time= double, typename Numeric=Time, int Dim=3, bool Safe=false
, typename Point= Eigen::Matrix<Numeric, Dim, 1> >
struct  curve_abc : std::unary_function<Time, Point>
{
	typedef Point 	point_t;
	typedef Time 	time_t;

/* Constructors - destructors */
	public:
	///\brief Constructor
	SPLINE_API curve_abc(){};

	///\brief Destructor
	SPLINE_API ~curve_abc(){};
/* Constructors - destructors */

/*Operations*/
	public:
	///  \brief Evaluation of the cubic spline at time t.
	///  \param t : the time when to evaluate the spine
	///  \param return : the value x(t)
	SPLINE_API virtual point_t operator()(time_t t) const = 0;
/*Operations*/

/*Helpers*/
	SPLINE_API virtual time_t MinBound() const = 0;
	SPLINE_API virtual time_t MaxBound() const = 0;
/*Helpers*/

	};
}
#endif //_STRUCT_CURVE_ABC

