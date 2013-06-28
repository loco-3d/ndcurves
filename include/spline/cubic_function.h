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

#include "curve_abc.h"

#include <stdexcept>

namespace spline
{
/// \class CubicFunction
/// \brief Represents a cubic spline defined on the interval 
/// [tBegin, tEnd]. It follows the equation
/// x(t) = a + b(t - t_min_) + c(t - t_min_)^2 + d(t - t_min_)^3 
///
template<typename Time= double, typename Numeric=Time, int Dim=3, bool Safe=false
, typename Point= Eigen::Matrix<Numeric, Dim, 1> >
struct cubic_function : public curve_abc<Time, Numeric, Dim, Safe, Point>
{
	typedef Point 	point_t;
  	typedef Time 	time_t;
  	typedef Numeric	num_t;
/* Constructors - destructors */
	public:
	///\brief Constructor
	cubic_function(point_t const& a, point_t const& b, point_t const& c, point_t const &d, time_t min, time_t max)
		:a_(a), b_(b), c_(c), d_(d), t_min_(min), t_max_(max)
	{
		if(t_min_ >= t_max_ && Safe)
		{
			std::out_of_range("TODO");
		}
	}

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
	public:
	///  \brief Evaluation of the cubic spline at time t.
	///  \param t : the time when to evaluate the spine
	///  \param return : the value x(t)
	virtual point_t operator()(time_t t) const 
	{
	    	if((t < t_min_ || t > t_max_) && Safe){ throw std::out_of_range("TODO");}
	    	time_t const dt (t-t_min_);
	    	return a_+ b_ * dt + c_ * dt*dt  + d_ * dt*dt*dt;
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
	const point_t a_, b_, c_ ,d_;
	const time_t t_min_, t_max_;
/*Attributes*/
	}; //class CubicFunction
}
#endif //_STRUCT_CUBICFUNCTION

