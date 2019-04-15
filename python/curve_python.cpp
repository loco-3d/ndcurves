#include "curve/bezier_curve.h"
#include "curve/polynom.h"
#include "curve/exact_cubic.h"
#include "curve/spline_deriv_constraint.h"
#include "curve/curve_constraint.h"
#include "curve/bezier_polynom_conversion.h"
#include "curve/bernstein.h"


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

typedef curve::bezier_curve  <real, real, 3, true, point_t> bezier3_t;
typedef curve::bezier_curve  <real, real, 6, true, point6_t> bezier6_t;
typedef curve::polynom  <real, real, 3, true, point_t, t_point_t> polynom_t;
typedef curve::exact_cubic  <real, real, 3, true, point_t, t_point_t> exact_cubic_t;
typedef polynom_t::coeff_t coeff_t;
typedef std::pair<real, point_t> waypoint_t;
typedef std::vector<waypoint_t, Eigen::aligned_allocator<point_t> > t_waypoint_t;

typedef curve::Bern<double> bernstein_t;


typedef curve::spline_deriv_constraint  <real, real, 3, true, point_t, t_point_t> spline_deriv_constraint_t;
typedef curve::curve_constraints<point_t> curve_constraints_t;
typedef curve::curve_constraints<point6_t> curve_constraints6_t;
/*** TEMPLATE SPECIALIZATION FOR PYTHON ****/

EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(bernstein_t)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(bezier3_t)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(bezier6_t)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(polynom_t)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(exact_cubic_t)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(curve_constraints_t)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(spline_deriv_constraint_t)

namespace curve
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

template <typename Bezier, typename PointList, typename T_Point>
Bezier* wrapBezierConstructorTemplate(const PointList& array, const real ub =1.)
{
    T_Point asVector = vectorFromEigenArray<PointList, T_Point>(array);
    return new Bezier(asVector.begin(), asVector.end(), ub);
}

template <typename Bezier, typename PointList, typename T_Point, typename CurveConstraints>
Bezier* wrapBezierConstructorConstraintsTemplate(const PointList& array, const CurveConstraints& constraints, const real ub =1.)
{
    T_Point asVector = vectorFromEigenArray<PointList, T_Point>(array);
    return new Bezier(asVector.begin(), asVector.end(), constraints, ub);
}

/*3D constructors */
bezier3_t* wrapBezierConstructor(const point_list_t& array)
{
    return wrapBezierConstructorTemplate<bezier3_t, point_list_t, t_point_t>(array) ;
}
bezier3_t* wrapBezierConstructorBounds(const point_list_t& array, const real ub)
{
    return wrapBezierConstructorTemplate<bezier3_t, point_list_t, t_point_t>(array, ub) ;
}
bezier3_t* wrapBezierConstructorConstraints(const point_list_t& array, const curve_constraints_t& constraints)
{
    return wrapBezierConstructorConstraintsTemplate<bezier3_t, point_list_t, t_point_t, curve_constraints_t>(array, constraints) ;
}
bezier3_t* wrapBezierConstructorBoundsConstraints(const point_list_t& array, const curve_constraints_t& constraints, const real ub)
{
    return wrapBezierConstructorConstraintsTemplate<bezier3_t, point_list_t, t_point_t, curve_constraints_t>(array, constraints, ub) ;
}
/*END 3D constructors */
/*6D constructors */
bezier6_t* wrapBezierConstructor6(const point_list6_t& array)
{
    return wrapBezierConstructorTemplate<bezier6_t, point_list6_t, t_point6_t>(array) ;
}
bezier6_t* wrapBezierConstructorBounds6(const point_list6_t& array, const real ub)
{
    return wrapBezierConstructorTemplate<bezier6_t, point_list6_t, t_point6_t>(array, ub) ;
}
bezier6_t* wrapBezierConstructor6Constraints(const point_list6_t& array, const curve_constraints6_t& constraints)
{
    return wrapBezierConstructorConstraintsTemplate<bezier6_t, point_list6_t, t_point6_t, curve_constraints6_t>(array, constraints) ;
}
bezier6_t* wrapBezierConstructorBounds6Constraints(const point_list6_t& array, const curve_constraints6_t& constraints, const real ub)
{
    return wrapBezierConstructorConstraintsTemplate<bezier6_t, point_list6_t, t_point6_t, curve_constraints6_t>(array, constraints, ub) ;
}
/*END 6D constructors */

polynom_t* wrapSplineConstructor(const coeff_t& array)
{
    return new polynom_t(array, 0., 1.);
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


spline_deriv_constraint_t* wrapSplineDerivConstraint(const coeff_t& array, const time_waypoints_t& time_wp, const curve_constraints_t& constraints)
{
    t_waypoint_t wps = getWayPoints(array, time_wp);
    return new spline_deriv_constraint_t(wps.begin(), wps.end(),constraints);
}

spline_deriv_constraint_t* wrapSplineDerivConstraintNoConstraints(const coeff_t& array, const time_waypoints_t& time_wp)
{
    t_waypoint_t wps = getWayPoints(array, time_wp);
    return new spline_deriv_constraint_t(wps.begin(), wps.end());
}

point_t get_init_vel(const curve_constraints_t& c)
{
    return c.init_vel;
}

point_t get_init_acc(const curve_constraints_t& c)
{
    return c.init_acc;
}

point_t get_end_vel(const curve_constraints_t& c)
{
    return c.end_vel;
}

point_t get_end_acc(const curve_constraints_t& c)
{
    return c.end_acc;
}

void set_init_vel(curve_constraints_t& c, const point_t& val)
{
    c.init_vel = val;
}

void set_init_acc(curve_constraints_t& c, const point_t& val)
{
    c.init_acc = val;
}

void set_end_vel(curve_constraints_t& c, const point_t& val)
{
    c.end_vel = val;
}

void set_end_acc(curve_constraints_t& c, const point_t& val)
{
    c.end_acc = val;
}



BOOST_PYTHON_MODULE(curves)
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
            //.def("__init__", make_constructor(&wrapBezierConstructor6Constraints))
            //.def("__init__", make_constructor(&wrapBezierConstructorBounds6Constraints))
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
    class_<bezier3_t>
        ("bezier3", no_init)
            .def("__init__", make_constructor(&wrapBezierConstructor))
            .def("__init__", make_constructor(&wrapBezierConstructorBounds))
            .def("__init__", make_constructor(&wrapBezierConstructorConstraints))
            .def("__init__", make_constructor(&wrapBezierConstructorBoundsConstraints))
            .def("min", &bezier3_t::min)
            .def("max", &bezier3_t::max)
            .def("__call__", &bezier3_t::operator())
            .def("derivate", &bezier3_t::derivate)
            .def("compute_derivate", &bezier3_t::compute_derivate)
            .def("compute_primitive", &bezier3_t::compute_primitive)
            .def("waypoints", &wayPointsToList<bezier3_t,3>)
            .def_readonly("degree", &bezier3_t::degree_)
            .def_readonly("nbWaypoints", &bezier3_t::size_)
        ;
    /** END bezier curve**/


    /** BEGIN spline curve function**/
    class_<polynom_t>("polynom",  init<const polynom_t::coeff_t, const real, const real >())
            .def("__init__", make_constructor(&wrapSplineConstructor))
            .def("min", &polynom_t::min)
            .def("max", &polynom_t::max)
            .def("__call__", &polynom_t::operator())
            .def("derivate", &polynom_t::derivate)
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


    /** BEGIN curve constraints**/
    class_<curve_constraints_t>
        ("curve_constraints", init<>())
            .add_property("init_vel", &get_init_vel, &set_init_vel)
            .add_property("init_acc", &get_init_acc, &set_init_acc)
            .add_property("end_vel", &get_end_vel, &set_end_vel)
            .add_property("end_acc", &get_end_acc, &set_end_acc)
        ;
    /** END curve constraints**/


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

    /** BEGIN bernstein polynom**/
    class_<bernstein_t>
        ("bernstein", init<const unsigned int, const unsigned int>())
            .def("__call__", &bernstein_t::operator())
        ;
    /** END bernstein polynom**/

    /** BEGIN Bezier to polynom conversion**/
    def("from_bezier", from_bezier<bezier3_t,polynom_t>);
    /** END Bezier to polynom conversion**/


}

} // namespace curve
