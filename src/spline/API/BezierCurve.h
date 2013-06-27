/**
* \file BezierCurve.h
* \brief class allowing to create a Bezier curve of dimension 1 <= n <= 3.
* \author Steve T.
* \version 0.1
* \date 06/17/2013
*/


#ifndef _CLASS_BEZIERCURVE
#define _CLASS_BEZIERCURVE

#include "Curve_ABC.h"

#include "Exports.h"
#include "MathDefs.h"

#include <memory>
#include <vector>

namespace spline
{
	/// \class BezierCurve
	/// \brief Represents a curve
	///
	class BezierCurve : public Curve_ABC
	{
/* Constructors - destructors */
	public:
		///\brief Constructor
		///\param points: the points parametering the Bezier curve
		///\TODO : sor far size above 3 is ignored
		SPLINE_API BezierCurve(const T_Vector&  /*points*/, const Real minBound=0, const Real maxBound=1);

		///\brief Destructor
		SPLINE_API ~BezierCurve();

	private:
		BezierCurve(const BezierCurve&);
		BezierCurve& operator=(const BezierCurve&);
/* Constructors - destructors */

/*Operations*/
	public:
	///  \brief Evaluation of the cubic spline at time t.
	///  \param t : the time when to evaluate the spine
	///  \param result : a reference to the Point set to the x(t)
	///  \param return : true if evaluation is successful, false if t is out of range
	SPLINE_API virtual bool Evaluate(const Real /*t*/, Vector3& /*result*/) const;
/*Operations*/
	
SPLINE_API virtual Real MinBound() const;
SPLINE_API virtual Real MaxBound() const;
/*Helpers*/

	public:
		const int size_;
		const Real minBound_, maxBound_;
	
	private:
		const T_Vector  pts_;
	};
}
#endif //_CLASS_BEZIERCURVE

