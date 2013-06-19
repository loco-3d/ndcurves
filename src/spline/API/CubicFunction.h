/**
* \file ExactCubic.h
* \brief Definition of a cubic spline.
* \author Steve T.
* \version 0.1
* \date 06/17/2013
*
* This file contains definitions for the CubicFunction class.
* It allows the creation and evaluation of natural 3D
* smooth cubic splines
*/


#ifndef _CLASS_CUBICFUNCTIONIMP
#define _CLASS_CUBICFUNCTIONIMP

#include "Exports.h"
#include "MathDefs.h"

namespace spline
{
	class SplineVisitor;

	/// \class CubicFunction
	/// \brief Represents a cubic spline defined on the interval 
	/// [tBegin, tEnd]. It follows the equation
	/// x(t) = a + b(t - tBegin) + c(t - tBegin)^2 + d(t - tBegin)^3 
	///
	class CubicFunction
	{
/* Constructors - destructors */
	public:
		///\brief Constructor
		SPLINE_API CubicFunction(const Vector3& /*a*/, const Vector3& /*b*/, const Vector3& /*c*/, const Vector3& /*d*/, const Real /*tBegin*/, const Real /*tEnd*/);

		///\brief Destructor
		SPLINE_API ~CubicFunction();

	private:
		CubicFunction(const CubicFunction&);
		CubicFunction& operator=(const CubicFunction&);
/* Constructors - destructors */

/*Operations*/
	public:
		///  \brief Evaluation of the cubic spline at time t.
		///  \param t : the time when to evaluate the spine
		///  \param result : a reference to the Point set to the x(t)
		SPLINE_API bool Evaluate(const Real /*t*/, Vector3& /*result*/) const;
/*Operations*/

/*Helpers*/
	public:
///  \brief Given a timestep dt, returns a set of values for the exact spline
///	 separated by this timestep
///  \param visitor : an object called for each timestep in the spline interval.
///  \param dt : the timestep
///  \param result : a reference to the Point set to the x(t)
SPLINE_API void Accept(SplineVisitor& /*visitor*/, Real /*dt*/) const;
/*Helpers*/

/*Attributes*/
	private:
		const Vector3 a_, b_, c_ ,d_;
		const Real tBegin_, tEnd_;
/*Attributes*/
	}; //class CubicFunction
}
#endif //_CLASS_CUBICFUNCTIONIMP