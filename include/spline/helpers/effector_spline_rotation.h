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


#ifndef _CLASS_EFFECTOR_SPLINE_ROTATION
#define _CLASS_EFFECTOR_SPLINE_ROTATION

#include "spline/helpers/effector_spline.h"
#include <Eigen/Geometry>

namespace spline
{
namespace helpers
{
typedef Eigen::Matrix<Numeric, 4, 1> quat_t;
typedef Eigen::Ref<quat_t> quat_ref_t;
typedef const Eigen::Ref<const quat_t> quat_ref_const_t;
typedef Eigen::Matrix<Numeric, 7, 1> config_t;

/// \class effector_spline_rotation
/// \brief Represents a trajectory for and end effector
/// uses the method effector_spline to create a spline trajectory.
/// Additionally, handles the rotation of the effector as follows:
/// does not rotate during the take off and landing phase,
/// then uses a SLERP algorithm to interpolate the rotation in the quaternion
/// space.
class effector_spline_rotation
{
    /* Constructors - destructors */
    public:
    /// \brief Constructor.
    /// Given a set of waypoints, and the normal vector of the start and
    /// ending positions, automatically create the spline such that:
    /// + init and end velocities / accelerations are 0.
    /// + the effector lifts and lands exactly in the direction of the specified normals
    /// \param wayPointsBegin   : an iterator pointing to the first element of a waypoint container
    /// \param wayPointsEnd     : an iterator pointing to the end           of a waypoint container
    /// \param to_quat          : 4D vector, quaternion indicating rotation at take off(x, y, z, w)
    /// \param land_quat        : 4D vector, quaternion indicating rotation at landing (x, y, z, w)
    /// \param lift_normal      : normal to be followed by end effector at take-off
    /// \param land_normal      : normal to be followed by end effector at landing
    /// \param lift_offset      : length of the straight line along normal at take-off
    /// \param land_offset      : length of the straight line along normal at landing
    /// \param lift_offset_duration : time travelled along straight line at take-off
    /// \param land_offset_duration : time travelled along straight line at landing
    ///
    template<typename In>
    effector_spline_rotation(In wayPointsBegin, In wayPointsEnd,
            quat_ref_const_t& to_quat=quat_t(0,0,0,1), quat_ref_const_t& land_quat=quat_t(0,0,0,1),
            const Point& lift_normal=Eigen::Vector3d::UnitZ(), const Point& land_normal=Eigen::Vector3d::UnitZ(),
            const Numeric lift_offset=0.02, const Numeric land_offset=0.02,
            const Time lift_offset_duration=0.02, const Time land_offset_duration=0.02)
        : spline_(effector_spline(wayPointsBegin, wayPointsEnd, lift_normal, land_normal, lift_offset, land_offset, lift_offset_duration, land_offset_duration))
        , to_quat_(to_quat.data()), land_quat_(land_quat.data())
        , time_lift_offset_(spline_->min()+lift_offset_duration)
        , time_land_offset_(spline_->max()-land_offset_duration)
    {
        // NOTHING
    }

    ~effector_spline_rotation() {delete spline_;}
    /* Constructors - destructors */

    /*Helpers*/
    public:
    Numeric min() const{return spline_->min();}
    Numeric max() const{return spline_->max();}
    /*Helpers*/

    /*Operations*/
    public:
    ///  \brief Evaluation of the effector position and rotation at time t.
    ///  \param t : the time when to evaluate the spline
    ///  \param return : a 7D vector; The 3 first values are the 3D position, the 4 last the
    ///  quaternion describing the rotation
    ///
    config_t operator()(const Numeric t) const
    {
        config_t res;
        res.head<3>() = (*spline_)(t);
        res.tail<4>() = interpolate_quat(t);
        return res;
    }

    ///  \brief Interpolates between two quaternions
    ///  \param t : the time when to evaluate the spline
    ///  \param quat : quaternion updated as the interpolation result
    ///
    quat_t interpolate_quat(time_t t) const
    {
        if(t<=time_lift_offset_) return quat_t(to_quat_.coeffs().data());
        if(t>=time_land_offset_) return quat_t(land_quat_.coeffs().data());
        //normalize u
        Numeric u = (t - time_lift_offset_) /(time_land_offset_ - time_lift_offset_);
        return quat_t(to_quat_.slerp(u, land_quat_).coeffs().data());
    }
    /*Operations*/

    /*Attributes*/
    public:
    const exact_cubic_t* spline_;
    const Eigen::Quaterniond to_quat_;
    const Eigen::Quaterniond land_quat_;
    const double time_lift_offset_;
    const double time_land_offset_;
    /*Attributes*/
};

} // namespace helpers
} // namespace spline
#endif //_CLASS_EFFECTOR_SPLINE_ROTATION

