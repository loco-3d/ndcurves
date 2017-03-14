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
typedef Eigen::Matrix<double, 3, 1, 0, 3, 1> ret_point_t;
typedef Eigen::VectorXd time_waypoints_t;
typedef Eigen::Matrix<real, 3, Eigen::Dynamic> point_list_t;
typedef std::vector<point_t,Eigen::aligned_allocator<point_t> >  t_point_t;
typedef std::pair<real, point_t> Waypoint;
typedef std::vector<Waypoint> T_Waypoint;

typedef spline::bezier_curve  <real, real, 3, true, point_t> bezier_t;
typedef spline::spline_curve  <real, real, 3, true, point_t, t_point_t> spline_curve_t;
typedef spline::exact_cubic  <real, real, 3, true, point_t, t_point_t> exact_cubic_t;
typedef spline_curve_t::coeff_t coeff_t;
typedef std::pair<real, point_t> waypoint_t;
typedef std::vector<waypoint_t, Eigen::aligned_allocator<point_t> > t_waypoint_t;


typedef spline::spline_deriv_constraint  <real, real, 3, true, point_t, t_point_t> spline_deriv_constraint_t;
typedef spline_deriv_constraint_t::spline_constraints spline_constraints_t;
/*** TEMPLATE SPECIALIZATION FOR PYTHON ****/

EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(bezier_t)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(spline_curve_t)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(exact_cubic_t)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(spline_constraints_t)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(spline_deriv_constraint_t)

namespace spline
{
using namespace boost::python;

t_point_t vectorFromEigenArray(const point_list_t& array)
{
    t_point_t res;
    for(int i =0;i<array.cols();++i)
        res.push_back(array.col(i));
    return res;
}

bezier_t* wrapBezierConstructor(const point_list_t& array)
{
    t_point_t asVector = vectorFromEigenArray(array);
    return new bezier_t(asVector.begin(), asVector.end());
}


bezier_t* wrapBezierConstructorBounds(const point_list_t& array, const real lb, const real ub)
{
    t_point_t asVector = vectorFromEigenArray(array);
    return new bezier_t(asVector.begin(), asVector.end(), lb, ub);
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
    eigenpy::enableEigenPySpecific<coeff_t,coeff_t>();
    /*eigenpy::exposeAngleAxis();
    eigenpy::exposeQuaternion();*/
    /** END eigenpy init**/

    /** BEGIN bezier curve**/
    class_<bezier_t>
        ("bezier", no_init)
            .def("__init__", make_constructor(&wrapBezierConstructor))
            .def("__init__", make_constructor(&wrapBezierConstructorBounds))
            .def("min", &bezier_t::min)
            .def("max", &bezier_t::max)
            .def("__call__", &bezier_t::operator())
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
