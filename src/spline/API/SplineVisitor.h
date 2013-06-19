/**
* \file SplineVisitor.h
* \brief Visitor for the Spline classes
* \author Steve T.
* \version 0.1
* \date 06/17/2013
*/


#ifndef _CLASS_SPLINE_VISITOR
#define _CLASS_SPLINE_VISITOR

namespace spline
{
	/// \class SplineVisitor
	/// \brief Represents a generic visitor for any Spline in the library.
	///
	class SPLINE_API SplineVisitor
	{
/* Constructors - destructors */
	public:
		///\brief Constructor
		SplineVisitor(){}

		///\brief Destructor
		~SplineVisitor(){}

	private:
		SplineVisitor(const SplineVisitor&);
		SplineVisitor& operator=(const SplineVisitor&);
/* Constructors - destructors */

/*Operations*/
	public:
	///  \brief Evaluation of the cubic spline at time t.
	///  \param time : the time associated to the given value
	///  \param value : the value x(time)
	virtual void Visit(const Real /*time*/, const Vector3& /*value*/) = 0;
/*Operations*/
	};
}
#endif //_CLASS_SPLINE_VISITOR