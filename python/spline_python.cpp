#include "spline/bezier_curve.h"
#include "spline/spline_curve.h"
#include "spline/exact_cubic.h"
#include "spline/spline_deriv_constraint.h"

#include <vector>

#include <eigenpy/memory.hpp>
#include <eigenpy/eigenpy.hpp>

#include <boost/python.hpp>

/*** TEMPLATE SPECIALIZATION FOR PYTHON ****/
typedef double real;
typedef Eigen::Vector3d point_t;
typedef Eigen::Matrix<double, 6, 1, 0, 6, 1> point6_t;
typedef Eigen::Matrix<double, 3, 1, 0, 3, 1> ret_point_t;
typedef Eigen::Matrix<double, 6, 1, 0, 6, 1> ret_point6_t;
typedef Eigen::VectorXd time_waypoints_t;
typedef Eigen::Matrix<real, 3, Eigen::Dynamic> point_list_t;
typedef Eigen::Matrix<real, 6, Eigen::Dynamic> point_list6_t;
typedef std::vector<point_t,Eigen::aligned_allocator<point_t> >  t_point_t;
typedef std::vector<point6_t,Eigen::aligned_allocator<point6_t> >  t_point6_t;
typedef std::pair<real, point_t> Waypoint;
typedef std::vector<Waypoint> T_Waypoint;
typedef std::pair<real, point6_t> Waypoint6;
typedef std::vector<Waypoint6> T_Waypoint6;

typedef spline::bezier_curve  <real, real, 3, true, point_t> bezier_t;
typedef spline::bezier_curve  <real, real, 6, true, point6_t> bezier6_t;
typedef spline::spline_curve  <real, real, 3, true, point_t, t_point_t> spline_curve_t;
typedef spline::exact_cubic  <real, real, 3, true, point_t, t_point_t> exact_cubic_t;
typedef spline_curve_t::coeff_t coeff_t;
typedef std::pair<real, point_t> waypoint_t;
typedef std::vector<waypoint_t, Eigen::aligned_allocator<point_t> > t_waypoint_t;


typedef spline::spline_deriv_constraint  <real, real, 3, true, point_t, t_point_t> spline_deriv_constraint_t;
typedef spline_deriv_constraint_t::spline_constraints spline_constraints_t;
/*** TEMPLATE SPECIALIZATION FOR PYTHON ****/

EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(bezier_t)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(bezier6_t)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(spline_curve_t)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(exact_cubic_t)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(spline_constraints_t)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(spline_deriv_constraint_t)

namespace spline
{
using namespace boost::python;
template <typename PointList, typename T_Point>
T_Point vectorFromEigenArray(const PointList& array)
{
    T_Point res;
    for(int i =0;i<array.cols();++i)
        res.push_back(array.col(i));
    return res;
}

bezier_t* wrapBezierConstructor(const point_list_t& array)
{
    t_point_t asVector = vectorFromEigenArray<point_list_t, t_point_t>(array);
    return new bezier_t(asVector.begin(), asVector.end());
}


bezier_t* wrapBezierConstructorBounds(const point_list_t& array, const real lb, const real ub)
{
    t_point_t asVector = vectorFromEigenArray<point_list_t, t_point_t>(array);
    return new bezier_t(asVector.begin(), asVector.end(), lb, ub);
}

bezier6_t* wrapBezierConstructor6(const point_list6_t& array)
{
    t_point6_t asVector = vectorFromEigenArray<point_list6_t, t_point6_t>(array);
    return new bezier6_t(asVector.begin(), asVector.end());
}


bezier6_t* wrapBezierConstructorBounds6(const point_list6_t& array, const real lb, const real ub)
{
    t_point6_t asVector = vectorFromEigenArray<point_list6_t, t_point6_t>(array);
    return new bezier6_t(asVector.begin(), asVector.end(), lb, ub);
}

spline_curve_t* wrapSplineConstructor(const coeff_t& array)
{
    return new spline_curve_t(array, 0., 1.);
}


t_waypoint_t getWayPoints(const coeff_t& array, const time_waypoints_t& time_wp)
{
    t_waypoint_t res;
    for(int i =0;i<array.cols();++i)
        res.push_back(std::make_pair(time_wp(i), array.col(i)));
    return res;
}

template <typename BezierType, int dim>
Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic> wayPointsToList(const BezierType& self)
{
    typedef typename BezierType::t_point_t t_point;
    typedef typename BezierType::t_point_t::const_iterator cit_point;
    const t_point& wps = self.waypoints();
    Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic> res (dim, wps.size());
    int col = 0;
    for(cit_point cit = wps.begin(); cit != wps.end(); ++cit, ++col)
        res.block<dim,1>(0,col) = *cit;
    return res;
}

exact_cubic_t* wrapExactCubicConstructor(const coeff_t& array, const time_waypoints_t& time_wp)
{
    t_waypoint_t wps = getWayPoints(array, time_wp);
    return new exact_cubic_t(wps.begin(), wps.end());
}


spline_deriv_constraint_t* wrapSplineDerivConstraint(const coeff_t& array, const time_waypoints_t& time_wp, const spline_constraints_t& constraints)
{
    t_waypoint_t wps = getWayPoints(array, time_wp);
    return new spline_deriv_constraint_t(wps.begin(), wps.end(),constraints);
}

spline_deriv_constraint_t* wrapSplineDerivConstraintNoConstraints(const coeff_t& array, const time_waypoints_t& time_wp)
{
    t_waypoint_t wps = getWayPoints(array, time_wp);
    return new spline_deriv_constraint_t(wps.begin(), wps.end());
}

point_t get_init_vel(const spline_constraints_t& c)
{
    return c.init_vel;
}

point_t get_init_acc(const spline_constraints_t& c)
{
    return c.init_acc;
}

point_t get_end_vel(const spline_constraints_t& c)
{
    return c.end_vel;
}

point_t get_end_acc(const spline_constraints_t& c)
{
    return c.end_acc;
}

void set_init_vel(spline_constraints_t& c, const point_t& val)
{
    c.init_vel = val;
}

void set_init_acc(spline_constraints_t& c, const point_t& val)
{
    c.init_acc = val;
}

void set_end_vel(spline_constraints_t& c, const point_t& val)
{
    c.end_vel = val;
}

void set_end_acc(spline_constraints_t& c, const point_t& val)
{
    c.end_acc = val;
}



BOOST_PYTHON_MODULE(spline)
{
    /** BEGIN eigenpy init**/
    eigenpy::enableEigenPy();

    eigenpy::enableEigenPySpecific<point_t,point_t>();
    eigenpy::enableEigenPySpecific<ret_point_t,ret_point_t>();
    eigenpy::enableEigenPySpecific<point_list_t,point_list_t>();
    eigenpy::enableEigenPySpecific<point6_t,point6_t>();
    eigenpy::enableEigenPySpecific<ret_point6_t,ret_point6_t>();
    eigenpy::enableEigenPySpecific<point_list6_t,point_list6_t>();
    eigenpy::enableEigenPySpecific<coeff_t,coeff_t>();
    /*eigenpy::exposeAngleAxis();
    eigenpy::exposeQuaternion();*/
    /** END eigenpy init**/

    /** BEGIN bezier curve 6**/
    class_<bezier6_t>
        ("bezier6", no_init)
            .def("__init__", make_constructor(&wrapBezierConstructor6))
            .def("__init__", make_constructor(&wrapBezierConstructorBounds6))
            .def("min", &bezier6_t::min)
            .def("max", &bezier6_t::max)
            .def("__call__", &bezier6_t::operator())
            .def("derivate", &bezier6_t::derivate)
            .def("compute_derivate", &bezier6_t::compute_derivate)
            .def("compute_primitive", &bezier6_t::compute_primitive)
            .def("waypoints", &wayPointsToList<bezier6_t,6>)
            .def_readonly("degree", &bezier6_t::degree_)
            .def_readonly("nbWaypoints", &bezier6_t::size_)
        ;
    /** END bezier curve**/

    /** BEGIN bezier curve**/
    class_<bezier_t>
        ("bezier", no_init)
            .def("__init__", make_constructor(&wrapBezierConstructor))
            .def("__init__", make_constructor(&wrapBezierConstructorBounds))
            .def("min", &bezier_t::min)
            .def("max", &bezier_t::max)
            .def("__call__", &bezier_t::operator())
            .def("derivate", &bezier_t::derivate)
            .def("compute_derivate", &bezier_t::compute_derivate)
            .def("compute_primitive", &bezier_t::compute_primitive)
            .def("waypoints", &wayPointsToList<bezier_t,3>)
            .def_readonly("degree", &bezier_t::degree_)
            .def_readonly("nbWaypoints", &bezier_t::size_)
        ;
    /** END bezier curve**/


    /** BEGIN spline curve function**/
    class_<spline_curve_t>("spline",  init<const spline_curve_t::coeff_t, const real, const real >())
            .def("__init__", make_constructor(&wrapSplineConstructor))
            .def("min", &spline_curve_t::min)
            .def("max", &spline_curve_t::max)
            .def("__call__", &spline_curve_t::operator())
            .def("derivate", &spline_curve_t::derivate)
        ;
    /** END cubic function**/


    /** BEGIN exact_cubic curve**/
    class_<exact_cubic_t>
        ("exact_cubic", no_init)
            .def("__init__", make_constructor(&wrapExactCubicConstructor))
            .def("min", &exact_cubic_t::min)
            .def("max", &exact_cubic_t::max)
            .def("__call__", &exact_cubic_t::operator())
            .def("derivate", &exact_cubic_t::derivate)
        ;
    /** END bezier curve**/


    /** BEGIN spline constraints**/
    class_<spline_constraints_t>
        ("spline_constraints", init<>())
            .add_property("init_vel", &get_init_vel, &set_init_vel)
            .add_property("init_acc", &get_init_acc, &set_init_acc)
            .add_property("end_vel", &get_end_vel, &set_end_vel)
            .add_property("end_acc", &get_end_acc, &set_end_acc)
        ;
    /** END spline constraints**/


    /** BEGIN spline_deriv_constraints**/
    class_<spline_deriv_constraint_t>
        ("spline_deriv_constraint", no_init)
            .def("__init__", make_constructor(&wrapSplineDerivConstraint))
            .def("__init__", make_constructor(&wrapSplineDerivConstraintNoConstraints))
            .def("min", &exact_cubic_t::min)
            .def("max", &exact_cubic_t::max)
            .def("__call__", &exact_cubic_t::operator())
            .def("derivate", &exact_cubic_t::derivate)
        ;
    /** END spline_deriv_constraints**/


}

} // namespace spline
