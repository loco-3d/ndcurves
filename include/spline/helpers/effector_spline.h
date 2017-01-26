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


#ifndef _CLASS_EFFECTORSPLINE
#define _CLASS_EFFECTORSPLINE

#include "spline/spline_deriv_constraint.h"

namespace spline
{
namespace helpers
{
typedef double Numeric;
typedef double Time;
typedef Eigen::Matrix<Numeric, 3, 1> Point;
typedef std::vector<Point,Eigen::aligned_allocator<Point> > T_Point;
typedef std::pair<double, Point> Waypoint;
typedef std::vector<Waypoint> T_Waypoint;
typedef spline_deriv_constraint<Time, Numeric, 3, true, Point, T_Point> spline_deriv_constraint_t;
typedef spline_deriv_constraint_t::spline_constraints spline_constraints_t;
typedef spline_deriv_constraint_t::exact_cubic_t exact_cubic_t;
typedef spline_deriv_constraint_t::t_spline_t t_spline_t;
typedef spline_deriv_constraint_t::spline_t spline_t;

Waypoint compute_offset(const Waypoint& source, const Point& normal, const Numeric offset, const Time time_offset)
{
    //compute time such that the equation from source to offsetpoint is necessarily a line
    Numeric norm = normal.norm();
    assert(norm>0.);
    return std::make_pair(source.first + time_offset,(source.second + normal / norm* offset));
}


spline_t make_end_spline(const Point& normal, const Point& from, const Numeric offset,
                         const Time init_time, const Time time_offset)
{
    // compute spline from land way point to end point
    // constraints are null velocity and acceleration
    Numeric norm = normal.norm();
    assert(norm>0.);
    Point n = normal / norm;
    Point d = offset / (time_offset*time_offset*time_offset) * -n;
    Point c = -3 * d * time_offset;
    Point b = -c * time_offset;

    T_Point points;
    points.push_back(from);
    points.push_back(b);
    points.push_back(c);
    points.push_back(d);

    return spline_t(points.begin(), points.end(), init_time, init_time+time_offset);
}

spline_constraints_t compute_required_offset_velocity_acceleration(const spline_t& end_spline, const Time time_offset)
{
    //computing end velocity: along landing normal and respecting time
    spline_constraints_t constraints;
    constraints.end_acc = end_spline.derivate(end_spline.min(),2);
    constraints.end_vel = end_spline.derivate(end_spline.min(),1);
    return constraints;
}


/// \brief Helper method to create a spline typically used to
/// guide the 3d trajectory of a robot end effector.
/// Given a set of waypoints, and the normal vector of the start and
/// ending positions, automatically create the spline such that:
/// + init and end velocities / accelerations are 0.
/// + the effector lifts and lands exactly in the direction of the specified normals
/// \param wayPointsBegin   : an iterator pointing to the first element of a waypoint container
/// \param wayPointsEnd     : an iterator pointing to the end           of a waypoint container
/// \param lift_normal      : normal to be followed by end effector at take-off
/// \param land_normal      : normal to be followed by end effector at landing
/// \param lift_offset      : length of the straight line along normal at take-off
/// \param land_offset      : length of the straight line along normal at landing
/// \param lift_offset_duration : time travelled along straight line at take-off
/// \param land_offset_duration : time travelled along straight line at landing
///
template<typename In>
exact_cubic_t* effector_spline(
        In wayPointsBegin, In wayPointsEnd, const Point& lift_normal=Eigen::Vector3d::UnitZ(), const Point& land_normal=Eigen::Vector3d::UnitZ(),
        const Numeric lift_offset=0.02, const Numeric land_offset=0.02,
        const Time lift_offset_duration=0.02, const Time land_offset_duration=0.02)
{
    T_Waypoint waypoints;
    const Waypoint& inPoint=*wayPointsBegin, endPoint =*(wayPointsEnd-1);
    waypoints.push_back(inPoint);
    //adding initial offset
    waypoints.push_back(compute_offset(inPoint, lift_normal ,lift_offset, lift_offset_duration));
    //inserting all waypoints but last
    waypoints.insert(waypoints.end(),wayPointsBegin+1,wayPointsEnd-1);
    //inserting waypoint to start landing
    const Waypoint& landWaypoint=compute_offset(endPoint, land_normal ,land_offset, -land_offset_duration);
    waypoints.push_back(landWaypoint);
    //specifying end velocity constraint such that landing will be in straight line
    spline_t end_spline=make_end_spline(land_normal,landWaypoint.second,land_offset,landWaypoint.first,land_offset_duration);
    spline_constraints_t constraints = compute_required_offset_velocity_acceleration(end_spline,land_offset_duration);
    spline_deriv_constraint_t all_but_end(waypoints.begin(), waypoints.end(),constraints);
    t_spline_t splines = all_but_end.subSplines_;
    splines.push_back(end_spline);
    return new exact_cubic_t(splines);
}
}
}
#endif //_CLASS_EFFECTORSPLINE

