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
#include "curve_constraint.h"

#include "MathDefs.h"

#include <functional>
#include <vector>

namespace spline
{
/// \class spline_deriv_constraint.
/// \brief Represents a set of cubic splines defining a continuous function 
/// crossing each of the waypoint given in its initialization. Additional constraints
/// are used to increase the order of the last spline, to start and finish
/// trajectory with user defined velocity and acceleration.
///
///
template<typename Time= double, typename Numeric=Time, std::size_t Dim=3, bool Safe=false,
         typename Point= Eigen::Matrix<Numeric, Dim, 1>,
         typename T_Point =std::vector<Point,Eigen::aligned_allocator<Point> >,
         typename SplineBase=polynom<Time, Numeric, Dim, Safe, Point, T_Point> >
struct spline_deriv_constraint : public exact_cubic<Time, Numeric, Dim, Safe, Point, T_Point, SplineBase>
{
    typedef Point 	point_t;
    typedef T_Point t_point_t;
    typedef Eigen::Matrix<Numeric, Eigen::Dynamic, Eigen::Dynamic> MatrixX;
    typedef Eigen::Matrix<Numeric, 3, 3> Matrix3;
    typedef Time 	time_t;
    typedef Numeric	num_t;
    typedef polynom<time_t, Numeric, Dim, Safe, point_t, t_point_t> spline_t;
    typedef exact_cubic<time_t, Numeric, Dim, Safe, point_t, t_point_t> exact_cubic_t;
    typedef typename std::vector<spline_t> t_spline_t;
    typedef typename t_spline_t::iterator it_spline_t;
    typedef typename t_spline_t::const_iterator cit_spline_t;
    typedef curve_constraints<point_t> spline_constraints;

	/* Constructors - destructors */
	public:
	///\brief Constructor
	///\param wayPointsBegin : an iterator pointing to the first element of a waypoint container
    ///\param wayPointsEnd   : an iterator pointing to the end           of a waypoint container
    ///\param constraints    : constraints on the init and end velocity / accelerations of the spline
    template<typename In>
    spline_deriv_constraint(In wayPointsBegin, In wayPointsEnd, const spline_constraints& constraints = spline_constraints())
        : exact_cubic_t(computeWayPoints<In>(wayPointsBegin, wayPointsEnd, constraints)) {}

	///\brief Destructor
    ~spline_deriv_constraint(){}

    ///\brief Copy Constructor
    spline_deriv_constraint(const spline_deriv_constraint& other)
        : exact_cubic_t(other.subSplines_) {}

    private:
    template<typename In>
    void compute_one_spline(In wayPointsBegin, In wayPointsNext, spline_constraints& constraints, t_spline_t& subSplines) const
    {
        const point_t& a0 = wayPointsBegin->second, a1 = wayPointsNext->second;
        const point_t& b0 = constraints.init_vel , c0 = constraints.init_acc / 2.;
        const num_t& init_t = wayPointsBegin->first, end_t = wayPointsNext->first;
        const num_t dt = end_t - init_t, dt_2 = dt * dt, dt_3 = dt_2 * dt;
        const point_t d0 = (a1 - a0 - b0 * dt - c0 * dt_2) / dt_3;
        subSplines.push_back(create_cubic<Time,Numeric,Dim,Safe,Point,T_Point>
                             (a0,b0,c0,d0,init_t, end_t));
        constraints.init_vel = subSplines.back().derivate(end_t, 1);
        constraints.init_acc = subSplines.back().derivate(end_t, 2);
    }

    template<typename In>
    void compute_end_spline(In wayPointsBegin, In wayPointsNext, spline_constraints& constraints, t_spline_t& subSplines) const
    {
        const point_t& a0 = wayPointsBegin->second, a1 = wayPointsNext->second;
        const point_t& b0 = constraints.init_vel, b1 = constraints.end_vel,
                       c0 = constraints.init_acc / 2., c1 = constraints.end_acc;
        const num_t& init_t = wayPointsBegin->first, end_t = wayPointsNext->first;
        const num_t dt = end_t - init_t, dt_2 = dt * dt, dt_3 = dt_2 * dt, dt_4 = dt_3 * dt, dt_5 = dt_4 * dt;
        //solving a system of four linear eq with 4 unknows: d0, e0
        const point_t alpha_0 = a1 - a0 - b0 *dt -    c0 * dt_2;
        const point_t alpha_1 = b1 -      b0     - 2 *c0 * dt;
        const point_t alpha_2 = c1 -               2 *c0;
        const num_t x_d_0 = dt_3, x_d_1 = 3*dt_2, x_d_2 = 6*dt;
        const num_t x_e_0 = dt_4, x_e_1 = 4*dt_3, x_e_2 = 12*dt_2;
        const num_t x_f_0 = dt_5, x_f_1 = 5*dt_4, x_f_2 = 20*dt_3;

        point_t d, e, f;
        MatrixX rhs = MatrixX::Zero(3,Dim);
        rhs.row(0) = alpha_0; rhs.row(1) = alpha_1; rhs.row(2) = alpha_2;
        Matrix3 eq  = Matrix3::Zero(3,3);
        eq(0,0) = x_d_0; eq(0,1) = x_e_0; eq(0,2) = x_f_0;
        eq(1,0) = x_d_1; eq(1,1) = x_e_1; eq(1,2) = x_f_1;
        eq(2,0) = x_d_2; eq(2,1) = x_e_2; eq(2,2) = x_f_2;
        rhs = eq.inverse().eval() * rhs;
        d = rhs.row(0); e = rhs.row(1); f = rhs.row(2);

        subSplines.push_back(create_quintic<Time,Numeric,Dim,Safe,Point,T_Point>
                             (a0,b0,c0,d,e,f, init_t, end_t));
    }

    private:
    template<typename In>
    t_spline_t computeWayPoints(In wayPointsBegin, In wayPointsEnd, const spline_constraints& constraints) const
    {
        std::size_t const size(std::distance(wayPointsBegin, wayPointsEnd));
        if(Safe && size < 1) throw; // TODO
        t_spline_t subSplines; subSplines.reserve(size-1);
        spline_constraints cons = constraints;
        In it(wayPointsBegin), next(wayPointsBegin), end(wayPointsEnd-1);
        ++next;
        for(std::size_t i(0); next != end; ++next, ++it, ++i)
            compute_one_spline<In>(it, next, cons, subSplines);
        compute_end_spline<In>(it, next,cons, subSplines);
        return subSplines;
    }

    private:
    /* Constructors - destructors */
};
}
#endif //_CLASS_CUBICZEROVELACC

