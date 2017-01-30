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
#include "spline/curve_abc.h"
#include <Eigen/Geometry>

namespace spline
{
namespace helpers
{

typedef Eigen::Matrix<Numeric, 4, 1> quat_t;
typedef Eigen::Ref<quat_t> quat_ref_t;
typedef const Eigen::Ref<const quat_t> quat_ref_const_t;
typedef Eigen::Matrix<Numeric, 7, 1> config_t;
typedef curve_abc<Time, Numeric, 4, false, quat_t> curve_abc_quat_t;
typedef std::pair<Numeric, quat_t > waypoint_quat_t;
typedef std::vector<waypoint_quat_t> t_waypoint_quat_t;
typedef Eigen::Matrix<Numeric, 1, 1> point_one_dim_t;
typedef spline_deriv_constraint <Numeric, Numeric, 1, false, point_one_dim_t> spline_deriv_constraint_one_dim;
typedef std::pair<Numeric, point_one_dim_t > waypoint_one_dim_t;
typedef std::vector<waypoint_one_dim_t> t_waypoint_one_dim_t;


class rotation_spline: public curve_abc_quat_t
{
    public:
     rotation_spline (quat_ref_const_t quat_from=quat_t(0,0,0,1), quat_ref_const_t quat_to=quat_t(0,0,0,1),
                               const double min = 0., const double max = 1.)
         : curve_abc_quat_t()
         , quat_from_(quat_from.data()), quat_to_(quat_to.data()), min_(min), max_(max)
         , time_reparam_(computeWayPoints())
     {}

    ~rotation_spline() {}

    /* Copy Constructors / operator=*/
    rotation_spline& operator=(const rotation_spline& from)
    {
        quat_from_ = from.quat_from_;
        quat_to_ = from.quat_to_;
        min_ =from.min_; max_ = from.max_;
        time_reparam_ = spline_deriv_constraint_one_dim(from.time_reparam_);
    }
    /* Copy Constructors / operator=*/

    quat_t operator()(const Numeric t) const
    {
        if(t<=min()) return quat_from_.coeffs();
        if(t>=max()) return quat_to_.coeffs();
        //normalize u
        Numeric u = (t - min()) /(max() - min());
        // reparametrize u
        return quat_from_.slerp(time_reparam_(u)[0], quat_to_).coeffs();
    }

    virtual quat_t derivate(time_t /*t*/, std::size_t /*order*/) const
    {
        throw std::runtime_error("TODO quaternion spline does not implement derivate");
    }

    spline_deriv_constraint_one_dim computeWayPoints() const
    {
        // initializing time reparametrization for spline
        t_waypoint_one_dim_t waypoints;
        waypoints.push_back(std::make_pair(0,point_one_dim_t::Zero()));
        waypoints.push_back(std::make_pair(1,point_one_dim_t::Ones()));
        return spline_deriv_constraint_one_dim(waypoints.begin(), waypoints.end());
    }

    virtual time_t min() const{return min_;}
    virtual time_t max() const{return max_;}

    public:
    Eigen::Quaterniond quat_from_; //const
    Eigen::Quaterniond quat_to_;   //const
    double min_;                   //const
    double max_;                   //const
    spline_deriv_constraint_one_dim time_reparam_; //const
};


typedef exact_cubic<Time, Numeric, 4, false, quat_t, std::vector<quat_t,Eigen::aligned_allocator<quat_t> >, rotation_spline> exact_cubic_quat_t;


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
        , quat_spline_(simple_quat_spline())
    {
        // NOTHING
    }

    /// \brief Constructor.
    /// Given a set of waypoints, and the normal vector of the start and
    /// ending positions, automatically create the spline such that:
    /// + init and end velocities / accelerations are 0.
    /// + the effector lifts and lands exactly in the direction of the specified normals
    /// \param wayPointsBegin       : an iterator pointing to the first element of a waypoint container
    /// \param wayPointsEnd         : an iterator pointing to the end           of a waypoint container
    /// \param quatWayPointsBegin   : en iterator pointing to the first element of a 4D vector (x, y, z, w) container of
    ///  quaternions indicating rotation at specific time steps.
    /// \param quatWayPointsEnd     : en iterator pointing to the end           of a 4D vector (x, y, z, w) container of
    ///  quaternions indicating rotation at specific time steps.
    /// \param lift_normal          : normal to be followed by end effector at take-off
    /// \param land_normal          : normal to be followed by end effector at landing
    /// \param lift_offset          : length of the straight line along normal at take-off
    /// \param land_offset          : length of the straight line along normal at landing
    /// \param lift_offset_duration : time travelled along straight line at take-off
    /// \param land_offset_duration : time travelled along straight line at landing
    ///
    template<typename In, typename InQuat>
    effector_spline_rotation(In wayPointsBegin, In wayPointsEnd,
            InQuat quatWayPointsBegin, InQuat quatWayPointsEnd,
            const Point& lift_normal=Eigen::Vector3d::UnitZ(), const Point& land_normal=Eigen::Vector3d::UnitZ(),
            const Numeric lift_offset=0.02, const Numeric land_offset=0.02,
            const Time lift_offset_duration=0.02, const Time land_offset_duration=0.02)
        : spline_(effector_spline(wayPointsBegin, wayPointsEnd, lift_normal, land_normal, lift_offset, land_offset, lift_offset_duration, land_offset_duration))
        , to_quat_((quatWayPointsBegin->second).data()), land_quat_(((quatWayPointsEnd-1)->second).data())
        , time_lift_offset_(spline_->min()+lift_offset_duration)
        , time_land_offset_(spline_->max()-land_offset_duration)
        , quat_spline_(quat_spline(quatWayPointsBegin, quatWayPointsEnd))
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

    public:
    quat_t interpolate_quat(const Numeric t) const
    {
        if(t<=time_lift_offset_) return to_quat_.coeffs();
        if(t>=time_land_offset_) return land_quat_.coeffs();
        return quat_spline_(t);
    }

    private:
    exact_cubic_quat_t simple_quat_spline() const
    {
        std::vector<rotation_spline> splines;
        splines.push_back(rotation_spline(to_quat_.coeffs(),land_quat_.coeffs(),time_lift_offset_, time_land_offset_));
        return exact_cubic_quat_t(splines);
    }

    template <typename InQuat>
    exact_cubic_quat_t quat_spline(InQuat quatWayPointsBegin, InQuat quatWayPointsEnd) const
    {
        if(std::distance(quatWayPointsBegin, quatWayPointsEnd) < 1)
            return simple_quat_spline();
        exact_cubic_quat_t::t_spline_t splines;
        InQuat it(quatWayPointsBegin);
        Time init = time_lift_offset_;
        Eigen::Quaterniond current_quat = to_quat_;
        for(; it != quatWayPointsEnd; ++it)
        {
            splines.push_back(rotation_spline(current_quat.coeffs(), it->second, init, it->first));
            current_quat = it->second;
            init = it->first;
        }
        splines.push_back(rotation_spline(current_quat.coeffs(), land_quat_.coeffs(), init, time_land_offset_));
        return exact_cubic_quat_t(splines);
    }

    /*Operations*/

    /*Attributes*/
    public:
    const exact_cubic_t* spline_;
    const Eigen::Quaterniond to_quat_;
    const Eigen::Quaterniond land_quat_;
    const double time_lift_offset_;
    const double time_land_offset_;
    const exact_cubic_quat_t quat_spline_;
    /*Attributes*/
};

} // namespace helpers
} // namespace spline
#endif //_CLASS_EFFECTOR_SPLINE_ROTATION

