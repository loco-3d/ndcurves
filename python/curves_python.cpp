#include "curves/bezier_curve.h"
#include "curves/polynomial.h"
#include "curves/exact_cubic.h"
#include "curves/curve_constraint.h"
#include "curves/curve_conversion.h"
#include "curves/bernstein.h"
#include "curves/cubic_hermite_spline.h"
#include "curves/piecewise_polynomial_curve.h"

#include "python_definitions.h"
#include "python_variables.h"

#include <vector>

#include <eigenpy/memory.hpp>
#include <eigenpy/eigenpy.hpp>

#include <boost/python.hpp>

/*** TEMPLATE SPECIALIZATION FOR PYTHON ****/
using namespace curves;
typedef double real;
typedef Eigen::Vector3d point_t;
typedef Eigen::Vector3d tangent_t;
typedef std::pair<point_t, tangent_t> pair_point_tangent_t;
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
typedef std::vector<pair_point_tangent_t,Eigen::aligned_allocator<pair_point_tangent_t> > t_pair_point_tangent_t;

typedef curves::bezier_curve  <real, real, 3, true, point_t> bezier3_t;
typedef curves::bezier_curve  <real, real, 6, true, point6_t> bezier6_t;
typedef curves::polynomial  <real, real, 3, true, point_t, t_point_t> polynomial_t;
typedef curves::piecewise_polynomial_curve <real, real, 3, true, point_t, t_point_t> piecewise_polynomial_curve_t;
typedef curves::exact_cubic  <real, real, 3, true, point_t, t_point_t> exact_cubic_t;
typedef curves::cubic_hermite_spline <real, real, 3, true, point_t> cubic_hermite_spline_t;
typedef polynomial_t::coeff_t coeff_t;
typedef std::pair<real, point_t> waypoint_t;
typedef std::vector<waypoint_t, Eigen::aligned_allocator<point_t> > t_waypoint_t;

typedef curves::Bern<double> bernstein_t;

typedef curves::curve_constraints<point_t> curve_constraints_t;
typedef curves::curve_constraints<point6_t> curve_constraints6_t;
/*** TEMPLATE SPECIALIZATION FOR PYTHON ****/

EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(bernstein_t)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(bezier3_t)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(bezier6_t)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(polynomial_t)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(exact_cubic_t)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(cubic_hermite_spline_t)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(curve_constraints_t)

namespace curves
{
using namespace boost::python;

/* Template constructor bezier */
template <typename Bezier, typename PointList, typename T_Point>
Bezier* wrapBezierConstructorTemplate(const PointList& array, const real T_min =0., const real T_max =1.)
{
    T_Point asVector = vectorFromEigenArray<PointList, T_Point>(array);
    return new Bezier(asVector.begin(), asVector.end(), T_min, T_max);
}

template <typename Bezier, typename PointList, typename T_Point, typename CurveConstraints>
Bezier* wrapBezierConstructorConstraintsTemplate(const PointList& array, const CurveConstraints& constraints, 
                                                 const real T_min =0., const real T_max =1.)
{
    T_Point asVector = vectorFromEigenArray<PointList, T_Point>(array);
    return new Bezier(asVector.begin(), asVector.end(), constraints, T_min, T_max);
}
/* End Template constructor bezier */

/*3D constructors bezier */
bezier3_t* wrapBezierConstructor3(const point_list_t& array)
{
    return wrapBezierConstructorTemplate<bezier3_t, point_list_t, t_point_t>(array) ;
}
bezier3_t* wrapBezierConstructorBounds3(const point_list_t& array, const real T_min, const real T_max)
{
    return wrapBezierConstructorTemplate<bezier3_t, point_list_t, t_point_t>(array, T_min, T_max) ;
}
bezier3_t* wrapBezierConstructor3Constraints(const point_list_t& array, const curve_constraints_t& constraints)
{
    return wrapBezierConstructorConstraintsTemplate<bezier3_t, point_list_t, t_point_t, curve_constraints_t>(array, constraints) ;
}
bezier3_t* wrapBezierConstructorBounds3Constraints(const point_list_t& array, const curve_constraints_t& constraints,
                                                   const real T_min, const real T_max)
{
    return wrapBezierConstructorConstraintsTemplate<bezier3_t, point_list_t, t_point_t, curve_constraints_t>(array, constraints, T_min, T_max) ;
}
/*END 3D constructors bezier */
/*6D constructors bezier */
bezier6_t* wrapBezierConstructor6(const point_list6_t& array)
{
    return wrapBezierConstructorTemplate<bezier6_t, point_list6_t, t_point6_t>(array) ;
}
bezier6_t* wrapBezierConstructorBounds6(const point_list6_t& array, const real T_min, const real T_max)
{
    return wrapBezierConstructorTemplate<bezier6_t, point_list6_t, t_point6_t>(array, T_min, T_max) ;
}
bezier6_t* wrapBezierConstructor6Constraints(const point_list6_t& array, const curve_constraints6_t& constraints)
{
    return wrapBezierConstructorConstraintsTemplate<bezier6_t, point_list6_t, t_point6_t, curve_constraints6_t>(array, constraints) ;
}
bezier6_t* wrapBezierConstructorBounds6Constraints(const point_list6_t& array, const curve_constraints6_t& constraints, const real T_min, const real T_max)
{
    return wrapBezierConstructorConstraintsTemplate<bezier6_t, point_list6_t, t_point6_t, curve_constraints6_t>(array, constraints, T_min, T_max) ;
}
/*END 6D constructors bezier */

/* Wrap Cubic hermite spline */
t_pair_point_tangent_t getPairsPointTangent(const point_list_t& points, const point_list_t& tangents)
{
    t_pair_point_tangent_t res;
    if (points.size() != tangents.size())
    {
        throw std::length_error("size of points and tangents must be the same");
    }
    for(int i =0;i<points.cols();++i) 
    {
        res.push_back(pair_point_tangent_t(tangents.col(i), tangents.col(i)));
    }
    return res;
}

cubic_hermite_spline_t* wrapCubicHermiteSplineConstructor(const point_list_t& points, const point_list_t& tangents, const time_waypoints_t& time_pts) 
{
    t_pair_point_tangent_t ppt = getPairsPointTangent(points, tangents);
    std::vector<real> time_control_pts;
    for( int i =0; i<time_pts.size(); ++i)
    {
        time_control_pts.push_back(time_pts[i]);
    }
    return new cubic_hermite_spline_t(ppt.begin(), ppt.end(), time_control_pts);
}
/* End wrap Cubic hermite spline */

/* Wrap polynomial */
polynomial_t* wrapSplineConstructor(const coeff_t& array)
{
    return new polynomial_t(array, 0., 1.);
}
/* End wrap polynomial */

/* Wrap piecewise polynomial curve */
piecewise_polynomial_curve_t* wrapPiecewisePolynomialCurveConstructor(const polynomial_t& pol)
{
    return new piecewise_polynomial_curve_t(pol);
}
/* end wrap piecewise polynomial curve */

/* Wrap exact cubic spline */
t_waypoint_t getWayPoints(const coeff_t& array, const time_waypoints_t& time_wp)
{
    t_waypoint_t res;
    for(int i =0;i<array.cols();++i) 
    {
        res.push_back(std::make_pair(time_wp(i), array.col(i)));
    }
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
    {
        res.block<dim,1>(0,col) = *cit;
    }
    return res;
}

exact_cubic_t* wrapExactCubicConstructor(const coeff_t& array, const time_waypoints_t& time_wp)
{
    t_waypoint_t wps = getWayPoints(array, time_wp);
    return new exact_cubic_t(wps.begin(), wps.end());
}

exact_cubic_t* wrapExactCubicConstructorConstraint(const coeff_t& array, const time_waypoints_t& time_wp, const curve_constraints_t& constraints)
{
    t_waypoint_t wps = getWayPoints(array, time_wp);
    return new exact_cubic_t(wps.begin(), wps.end(), constraints);
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
/* End wrap exact cubic spline */


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
            .def("__init__", make_constructor(&wrapBezierConstructor3))
            .def("__init__", make_constructor(&wrapBezierConstructorBounds3))
            .def("__init__", make_constructor(&wrapBezierConstructor3Constraints))
            .def("__init__", make_constructor(&wrapBezierConstructorBounds3Constraints))
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


    /** BEGIN variable points bezier curve**/
    class_<LinearControlPointsHolder>
        ("LinearWaypoint", no_init)
        .def_readonly("A", &LinearControlPointsHolder::A)
        .def_readonly("b", &LinearControlPointsHolder::b)
        ;

    class_<LinearBezierVector>
        ("bezierVarVector", no_init)
        .def_readonly("size", &LinearBezierVector::size)
        .def("at", &LinearBezierVector::at, return_value_policy<manage_new_object>())
        ;

    class_<bezier_linear_variable_t>
        ("bezierVar", no_init)
            .def("__init__", make_constructor(&wrapBezierLinearConstructor))
            .def("__init__", make_constructor(&wrapBezierLinearConstructorBounds))
            .def("min", &bezier_linear_variable_t::min)
            .def("max", &bezier_linear_variable_t::max)
            //.def("__call__", &bezier_linear_control_t::operator())
            .def("derivate", &bezier_linear_variable_t::derivate)
            .def("compute_derivate", &bezier_linear_variable_t::compute_derivate)
            .def("compute_primitive", &bezier_linear_variable_t::compute_primitive)
            .def("split", &split, return_value_policy<manage_new_object>())
            .def("waypoints", &wayPointsToLists, return_value_policy<manage_new_object>())
            .def_readonly("degree", &bezier_linear_variable_t::degree_)
            .def_readonly("nbWaypoints", &bezier_linear_variable_t::size_)
        ;
    /** END variable points bezier curve**/


    /** BEGIN polynomial curve function**/
    class_<polynomial_t>("polynomial",  init<const polynomial_t::coeff_t, const real, const real >())
            .def("__init__", make_constructor(&wrapSplineConstructor))
            .def("min", &polynomial_t::min)
            .def("max", &polynomial_t::max)
            .def("__call__", &polynomial_t::operator())
            .def("derivate", &polynomial_t::derivate)
        ;
    /** END polynomial function**/

    /** BEGIN piecewise polynomial curve function **/
    class_<piecewise_polynomial_curve_t>
        ("piecewise_polynomial_curve", no_init)
            .def("__init__", make_constructor(&wrapPiecewisePolynomialCurveConstructor))
            .def("min", &piecewise_polynomial_curve_t::min)
            .def("max", &piecewise_polynomial_curve_t::max)
            .def("__call__", &piecewise_polynomial_curve_t::operator())
            .def("derivate", &piecewise_polynomial_curve_t::derivate)
            .def("add_polynomial_curve", &piecewise_polynomial_curve_t::add_polynomial_curve)
            .def("is_continuous", &piecewise_polynomial_curve_t::is_continuous)
        ;
    /** END piecewise polynomial curve function **/


    /** BEGIN exact_cubic curve**/
    class_<exact_cubic_t>
        ("exact_cubic", no_init)
            .def("__init__", make_constructor(&wrapExactCubicConstructor))
            .def("__init__", make_constructor(&wrapExactCubicConstructorConstraint))
            .def("min", &exact_cubic_t::min)
            .def("max", &exact_cubic_t::max)
            .def("__call__", &exact_cubic_t::operator())
            .def("derivate", &exact_cubic_t::derivate)
            .def("getNumberSplines", &exact_cubic_t::getNumberSplines)
            .def("getSplineAt", &exact_cubic_t::getSplineAt)
        ;
    /** END exact_cubic curve**/


    /** BEGIN cubic_hermite_spline **/
    class_<cubic_hermite_spline_t>
        ("cubic_hermite_spline", no_init)
            .def("__init__", make_constructor(&wrapCubicHermiteSplineConstructor))
            .def("min", &cubic_hermite_spline_t::min)
            .def("max", &cubic_hermite_spline_t::max)
            .def("__call__", &cubic_hermite_spline_t::operator())
            .def("derivate", &cubic_hermite_spline_t::derivate)
        ;
    /** END cubic_hermite_spline **/


    /** BEGIN curve constraints**/
    class_<curve_constraints_t>
        ("curve_constraints", init<>())
            .add_property("init_vel", &get_init_vel, &set_init_vel)
            .add_property("init_acc", &get_init_acc, &set_init_acc)
            .add_property("end_vel", &get_end_vel, &set_end_vel)
            .add_property("end_acc", &get_end_acc, &set_end_acc)
        ;
    /** END curve constraints**/

    /** BEGIN bernstein polynomial**/
    class_<bernstein_t>
        ("bernstein", init<const unsigned int, const unsigned int>())
            .def("__call__", &bernstein_t::operator())
        ;
    /** END bernstein polynomial**/

    /** BEGIN curves conversion**/
    def("polynomial_from_bezier", polynomial_from_curve<polynomial_t,bezier3_t>);
    def("polynomial_from_hermite", polynomial_from_curve<polynomial_t,cubic_hermite_spline_t>);
    def("bezier_from_hermite", bezier_from_curve<bezier3_t,cubic_hermite_spline_t>);
    def("bezier_from_polynomial", bezier_from_curve<bezier3_t,polynomial_t>);
    def("hermite_from_bezier", hermite_from_curve<cubic_hermite_spline_t, bezier3_t>);
    def("hermite_from_polynomial", hermite_from_curve<cubic_hermite_spline_t, polynomial_t>);
    /** END curves conversion**/


}

} // namespace curves
