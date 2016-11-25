/**
* \file exact_cubic.h
* \brief class allowing to create an Exact cubic spline.
* \author Steve T.
* \version 0.1
* \date 06/17/2013
*
* This file contains definitions for the ExactCubic class.
* Given a set of waypoints (x_i*) and timestep (t_i), it provides the unique set of
* cubic splines fulfulling those 4 restrictions :
* - x_i(t_i) = x_i* ; this means that the curve passes trough each waypoint
* - x_i(t_i+1) = x_i+1* ;
* - its derivative is continous at t_i+1
* - its acceleration is continous at t_i+1
* more details in paper "Task-Space Trajectories via Cubic Spline Optimization"
* By J. Zico Kolter and Andrew Y.ng (ICRA 2009)
*/


#ifndef _CLASS_CUBICZEROVELACC
#define _CLASS_CUBICZEROVELACC

#include "exact_cubic.h"

#include "MathDefs.h"

#include <functional>
#include <vector>

namespace spline
{
/// \class cubic_zero_vel.
/// \brief Represents a set of cubic splines defining a continuous function 
/// crossing each of the waypoint given in its initialization. Additional constraints
/// are used to increase the order of the last and first splines, to start and finish
/// trajectory with zero velocity and acceleration. Thus the first and last splines
///
///
template<typename Time= double, typename Numeric=Time, std::size_t Dim=3, bool Safe=false,
         typename Point= Eigen::Matrix<Numeric, Dim, 1>,
         typename T_Point =std::vector<Point,Eigen::aligned_allocator<Point> > >
struct cubic_zero_vel : public exact_cubic<Time, Numeric, Dim, Safe, Point, T_Point>
{
    typedef Point 	point_t;
    typedef T_Point t_point_t;
    typedef Eigen::Matrix<Numeric, Eigen::Dynamic, Eigen::Dynamic> MatrixX;
    typedef Time 	time_t;
    typedef Numeric	num_t;
    typedef spline_curve<time_t, Numeric, Dim, Safe, point_t, t_point_t> spline_t;
    typedef exact_cubic<Time, Numeric, Dim, Safe, Point> exact_cubic_t;
    typedef typename std::vector<spline_t> t_spline_t;
    typedef typename t_spline_t::iterator it_spline_t;
    typedef typename t_spline_t::const_iterator cit_spline_t;

	/* Constructors - destructors */
	public:
	///\brief Constructor
	///\param wayPointsBegin : an iterator pointing to the first element of a waypoint container
	///\param wayPointsEns   : an iterator pointing to the end           of a waypoint container
	template<typename In>
    cubic_zero_vel(In wayPointsBegin, In wayPointsEnd)
        : exact_cubic_t(wayPointsBegin, wayPointsEnd)
	{
        // NOTHING
	}

	///\brief Destructor
    ~cubic_zero_vel(){}

	private:
    cubic_zero_vel(const cubic_zero_vel&);
    cubic_zero_vel& operator=(const cubic_zero_vel&);
    /* Constructors - destructors */
};
}
#endif //_CLASS_CUBICZEROVELACC

