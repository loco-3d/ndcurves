/**
* \file Curve_ABC.h
* \brief class allowing to create an Exact cubic spline.
* \author Steve T.
* \version 0.1
* \date 06/17/2013
*
* Interface for a curve
*/


#ifndef _CLASS_CURVEABC
#define _CLASS_CURVEABC

#include "SplineVisitor.h"

#include "Exports.h"
#include "MathDefs.h"

#include <memory>
#include <vector>

namespace spline
{
	/// \class Curve_ABC
	/// \brief Represents a curve
	///
	class Curve_ABC
	{
/* Constructors - destructors */
	public:
		///\brief Constructor
		SPLINE_API Curve_ABC(){};

		///\brief Destructor
		SPLINE_API ~Curve_ABC(){};

	private:
		Curve_ABC(const Curve_ABC&);
		Curve_ABC& operator=(const Curve_ABC&);
/* Constructors - destructors */

/*Operations*/
	public:
	///  \brief Evaluation of the cubic spline at time t.
	///  \param t : the time when to evaluate the spine
	///  \param result : a reference to the Point set to the x(t)
	///  \param return : true if evaluation is successful, false if t is out of range
	SPLINE_API virtual bool Evaluate(const Real /*t*/, Vector3& /*result*/) const = 0;
/*Operations*/

/*Helpers*/
	public:
///  \brief Given a timestep dt, returns a set of values for the exact spline
///	 separated by this timestep
///  \param visitor : an object called for each timestep in the spline interval.
///  \param dt : the timestep
///  \param result : a reference to the Point set to the x(t)
SPLINE_API void Accept(SplineVisitor& visitor, Real dt) const
{
	for(Real ti = MinBound(); ti <=  MaxBound(); ti = ti + dt)
	{
		Vector3 res; Evaluate(ti, res);
		visitor.Visit(ti, res);
	}
}

SPLINE_API virtual Real MinBound() const = 0;
SPLINE_API virtual Real MaxBound() const = 0;
/*Helpers*/

	};
}
#endif //_CLASS_EXACTCUBIC

