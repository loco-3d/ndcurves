#include "curves/bezier_curve.h"
#include "curves/polynomial.h"
#include "curves/exact_cubic.h"
#include "curves/curve_constraint.h"
#include "curves/curve_conversion.h"
#include "curves/bernstein.h"
#include "curves/cubic_hermite_spline.h"
#include "curves/piecewise_curve.h"

#include "python_definitions.h"
#include "python_variables.h"
#include "archive_python_binding.h"

#include <vector>

#include <eigenpy/memory.hpp>
#include <eigenpy/eigenpy.hpp>
#include <Eigen/Dense>

#include <boost/python.hpp>

/*** TEMPLATE SPECIALIZATION FOR PYTHON ****/
using namespace curves;
typedef double real;
typedef Eigen::VectorXd time_waypoints_t;
// 3D
typedef Eigen::Vector3d point3_t;
typedef Eigen::Matrix<double, 3, 1, 0, 3, 1> ret_point3_t;
typedef std::pair<point3_t, point3_t> pair_point3_tangent_t;
typedef Eigen::Matrix<real, 3, Eigen::Dynamic> point3_list_t;
typedef std::vector<point3_t,Eigen::aligned_allocator<point3_t> >  t_point3_t;
typedef std::vector<pair_point3_tangent_t,Eigen::aligned_allocator<pair_point3_tangent_t> > t_pair_point3_tangent_t;
// XD
typedef Eigen::VectorXd pointX_t;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1, 0, Eigen::Dynamic, 1> ret_pointX_t;
typedef std::pair<pointX_t, pointX_t> pair_pointX_tangent_t;
typedef Eigen::MatrixXd pointX_list_t;
typedef std::vector<pointX_t,Eigen::aligned_allocator<point3_t> >  t_pointX_t;
typedef std::vector<pair_pointX_tangent_t,Eigen::aligned_allocator<pair_pointX_tangent_t> > t_pair_pointX_tangent_t;
// Else
typedef std::pair<real, point3_t> Waypoint;
typedef std::vector<Waypoint> T_Waypoint;
typedef std::vector<Waypoint6> T_Waypoint6;

// Dynamic dim
typedef curves::cubic_hermite_spline <real, real, true, pointX_t> cubic_hermite_spline_t;

// 3D
typedef curves::bezier_curve  <real, real, 3, true, point3_t> bezier3_t;
typedef curves::polynomial  <real, real, 3, true, point3_t, t_point3_t> polynomial3_t;
typedef curves::exact_cubic  <real, real, 3, true, point3_t, t_point3_t> exact_cubic3_t;
typedef polynomial3_t::coeff_t coeff_t;
typedef std::pair<real, point3_t> waypoint3_t;
typedef std::vector<waypoint3_t, Eigen::aligned_allocator<point3_t> > t_waypoint3_t;

typedef curves::piecewise_curve <real, real, 3, true, point3_t, t_point3_t, polynomial3_t> piecewise_polynomial3_curve_t;
typedef curves::piecewise_curve <real, real, 3, true, point3_t, t_point3_t, bezier3_t> piecewise_bezier3_curve_t;
typedef curves::piecewise_curve <real, real, 3, true, point3_t, t_point3_t, cubic_hermite_spline_t> piecewise_cubic_hermite_curve_t;

typedef curves::Bern<double> bernstein_t;

typedef curves::curve_constraints<point3_t> curve_constraints3_t;
/*** TEMPLATE SPECIALIZATION FOR PYTHON ****/

EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(bernstein_t)

EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(bezier3_t)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(polynomial3_t)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(exact_cubic3_t)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(cubic_hermite_spline_t)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(curve_constraints3_t)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(piecewise_polynomial3_curve_t)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(piecewise_bezier3_curve_t)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(piecewise_cubic_hermite_curve_t)

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
  bezier3_t* wrapBezierConstructor3(const point3_list_t& array)
  {
    return wrapBezierConstructorTemplate<bezier3_t, point3_list_t, t_point3_t>(array) ;
  }
  bezier3_t* wrapBezierConstructorBounds3(const point3_list_t& array, const real T_min, const real T_max)
  {
    return wrapBezierConstructorTemplate<bezier3_t, point3_list_t, t_point3_t>(array, T_min, T_max) ;
  }
  bezier3_t* wrapBezierConstructor3Constraints(const point3_list_t& array, const curve_constraints3_t& constraints)
  {
    return wrapBezierConstructorConstraintsTemplate<bezier3_t, point3_list_t, t_point3_t, curve_constraints3_t>(array, constraints) ;
  }
  bezier3_t* wrapBezierConstructorBounds3Constraints(const point3_list_t& array, const curve_constraints3_t& constraints,
                             const real T_min, const real T_max)
  {
  return wrapBezierConstructorConstraintsTemplate<bezier3_t, point3_list_t, t_point3_t, curve_constraints3_t>(array, constraints, 
                                                                                                           T_min, T_max) ;
  }
  /*END 3D constructors bezier */

  /* Wrap Cubic hermite spline */
  t_pair_pointX_tangent_t getPairsPointTangent(const pointX_list_t& points, const pointX_list_t& tangents)
  {
    t_pair_pointX_tangent_t res;
    if (points.size() != tangents.size())
    {
      throw std::length_error("size of points and tangents must be the same");
    }
    for(int i =0;i<points.cols();++i)
    {
      res.push_back(pair_pointX_tangent_t(tangents.col(i), tangents.col(i)));
    }
    return res;
  }

  cubic_hermite_spline_t* wrapCubicHermiteSplineConstructor(const pointX_list_t& points, const pointX_list_t& tangents, 
                                                            const time_waypoints_t& time_pts)
  {
    t_pair_pointX_tangent_t ppt = getPairsPointTangent(points, tangents);
    std::vector<real> time_control_pts;
    for( int i =0; i<time_pts.size(); ++i)
    {
      time_control_pts.push_back(time_pts[i]);
    }
    return new cubic_hermite_spline_t(ppt.begin(), ppt.end(), time_control_pts);
  }
  /* End wrap Cubic hermite spline */

  /* Wrap polynomial */
  polynomial3_t* wrapPolynomial3Constructor1(const coeff_t& array, const real min, const real max)
  {
    return new polynomial3_t(array, min, max);
  }
  polynomial3_t* wrapPolynomial3Constructor2(const coeff_t& array)
  {
    return new polynomial3_t(array, 0., 1.);
  }
  /* End wrap polynomial */

  /* Wrap piecewise curve */
  piecewise_polynomial3_curve_t* wrapPiecewisePolynomial3CurveConstructor(const polynomial3_t& pol)
  {
    return new piecewise_polynomial3_curve_t(pol);
  }
  piecewise_polynomial3_curve_t* wrapPiecewisePolynomial3CurveEmptyConstructor()
  {
    return new piecewise_polynomial3_curve_t();
  }
  piecewise_bezier3_curve_t* wrapPiecewiseBezier3CurveConstructor(const bezier3_t& bc)
  {
    return new piecewise_bezier3_curve_t(bc);
  }
  piecewise_bezier3_curve_t* wrapPiecewiseBezier3CurveEmptyConstructor()
  {
    return new piecewise_bezier3_curve_t();
  }
  piecewise_cubic_hermite_curve_t* wrapPiecewiseCubicHermite3CurveConstructor(const cubic_hermite_spline_t& ch)
  {
    return new piecewise_cubic_hermite_curve_t(ch);
  }
  piecewise_cubic_hermite_curve_t* wrapPiecewiseCubicHermite3CurveEmptyConstructor()
  {
    return new piecewise_cubic_hermite_curve_t();
  }
  /* end wrap piecewise polynomial curve */

  /* Wrap exact cubic spline */
  t_waypoint3_t getWayPoints(const coeff_t& array, const time_waypoints_t& time_wp)
  {
    t_waypoint3_t res;
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

  exact_cubic3_t* wrapExactCubic3Constructor(const coeff_t& array, const time_waypoints_t& time_wp)
  {
    t_waypoint3_t wps = getWayPoints(array, time_wp);
    return new exact_cubic3_t(wps.begin(), wps.end());
  }

  exact_cubic3_t* wrapExactCubic3ConstructorConstraint(const coeff_t& array, const time_waypoints_t& time_wp, 
                                                     const curve_constraints3_t& constraints)
  {
    t_waypoint3_t wps = getWayPoints(array, time_wp);
    return new exact_cubic3_t(wps.begin(), wps.end(), constraints);
  }

  point3_t get_init_vel(const curve_constraints3_t& c)
  {
    return c.init_vel;
  }

  point3_t get_init_acc(const curve_constraints3_t& c)
  {
    return c.init_acc;
  }

  point3_t get_end_vel(const curve_constraints3_t& c)
  {
    return c.end_vel;
  }

  point3_t get_end_acc(const curve_constraints3_t& c)
  {
    return c.end_acc;
  }

  void set_init_vel(curve_constraints3_t& c, const point3_t& val)
  {
    c.init_vel = val;
  }

  void set_init_acc(curve_constraints3_t& c, const point3_t& val)
  {
    c.init_acc = val;
  }

  void set_end_vel(curve_constraints3_t& c, const point3_t& val)
  {
    c.end_vel = val;
  }

  void set_end_acc(curve_constraints3_t& c, const point3_t& val)
  {
    c.end_acc = val;
  }
  /* End wrap exact cubic spline */


  // TO DO : Replace all load and save function for serialization in class by using 
  //         SerializableVisitor in archive_python_binding.
  BOOST_PYTHON_MODULE(curves)
  {
    /** BEGIN eigenpy init**/
    eigenpy::enableEigenPy();
    eigenpy::enableEigenPySpecific<pointX_t,pointX_t>();
    eigenpy::enableEigenPySpecific<pointX_list_t,pointX_list_t>();
    eigenpy::enableEigenPySpecific<point3_t,point3_t>();
    eigenpy::enableEigenPySpecific<point3_list_t,point3_list_t>();
    eigenpy::enableEigenPySpecific<coeff_t,coeff_t>();
    /*eigenpy::exposeAngleAxis();
    eigenpy::exposeQuaternion();*/
    /** END eigenpy init**/
    /** BEGIN bezier curve**/
    class_<bezier3_t>("bezier3", init<>())
      .def("__init__", make_constructor(&wrapBezierConstructor3))
      .def("__init__", make_constructor(&wrapBezierConstructorBounds3))
      .def("__init__", make_constructor(&wrapBezierConstructor3Constraints))
      .def("__init__", make_constructor(&wrapBezierConstructorBounds3Constraints))
      .def("min", &bezier3_t::min)
      .def("max", &bezier3_t::max)
      .def("dim", &bezier3_t::dim)
      .def("__call__", &bezier3_t::operator())
      .def("derivate", &bezier3_t::derivate)
      .def("compute_derivate", &bezier3_t::compute_derivate)
      .def("compute_primitive", &bezier3_t::compute_primitive)
      .def("waypoints", &wayPointsToList<bezier3_t,3>)
      .def("saveAsText", &bezier3_t::saveAsText<bezier3_t>,bp::args("filename"),"Saves *this inside a text file.")
      .def("loadFromText",&bezier3_t::loadFromText<bezier3_t>,bp::args("filename"),"Loads *this from a text file.")
      .def("saveAsXML",&bezier3_t::saveAsXML<bezier3_t>,bp::args("filename","tag_name"),"Saves *this inside a XML file.")
      .def("loadFromXML",&bezier3_t::loadFromXML<bezier3_t>,bp::args("filename","tag_name"),"Loads *this from a XML file.")
      .def("saveAsBinary",&bezier3_t::saveAsBinary<bezier3_t>,bp::args("filename"),"Saves *this inside a binary file.")
      .def("loadFromBinary",&bezier3_t::loadFromBinary<bezier3_t>,bp::args("filename"),"Loads *this from a binary file.")
      .def_readonly("degree", &bezier3_t::degree_)
      .def_readonly("nbWaypoints", &bezier3_t::size_)
      //.def(SerializableVisitor<bezier3_t>())
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
    class_<bezier_linear_variable_t>("bezierVar", no_init)
      .def("__init__", make_constructor(&wrapBezierLinearConstructor))
      .def("__init__", make_constructor(&wrapBezierLinearConstructorBounds))
      .def("min", &bezier_linear_variable_t::min)
      .def("max", &bezier_linear_variable_t::max)
      .def("dim", &bezier_linear_variable_t::dim)
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
    class_<polynomial3_t>("polynomial3",  init<>())
      .def("__init__", make_constructor(&wrapPolynomial3Constructor1),
           "Create polynomial spline from an Eigen matrix of coefficient defined for t \in [min,max]."
           " The matrix should contain one coefficient per column, from the zero order coefficient,up to the highest order."
           " Spline order is given by the number of the columns -1.")
      .def("__init__", make_constructor(&wrapPolynomial3Constructor2),
           "Create polynomial spline from an Eigen matrix of coefficient defined for t \in [0,1]."
           " The matrix should contain one coefficient per column, from the zero order coefficient,up to the highest order."
           " Spline order is given by the number of the columns -1.")
      .def("min", &polynomial3_t::min, "Get the LOWER bound on interval definition of the curve.")
      .def("max", &polynomial3_t::max,"Get the HIGHER bound on interval definition of the curve.")
      .def("dim", &polynomial3_t::dim)
      .def("__call__", &polynomial3_t::operator(),"Evaluate the spline at the given time.")
      .def("derivate", &polynomial3_t::derivate,"Evaluate the derivative of order N of curve at time t.",args("self","t","N"))
      .def("compute_derivate", &polynomial3_t::compute_derivate,"Compute derivative of order N of curve at time t.")
      .def("saveAsText", &polynomial3_t::saveAsText<polynomial3_t>,bp::args("filename"),"Saves *this inside a text file.")
      .def("loadFromText",&polynomial3_t::loadFromText<polynomial3_t>,bp::args("filename"),"Loads *this from a text file.")
      .def("saveAsXML",&polynomial3_t::saveAsXML<polynomial3_t>,bp::args("filename","tag_name"),"Saves *this inside a XML file.")
      .def("loadFromXML",&polynomial3_t::loadFromXML<polynomial3_t>,bp::args("filename","tag_name"),"Loads *this from a XML file.")
      .def("saveAsBinary",&polynomial3_t::saveAsBinary<polynomial3_t>,bp::args("filename"),"Saves *this inside a binary file.")
      .def("loadFromBinary",&polynomial3_t::loadFromBinary<polynomial3_t>,bp::args("filename"),"Loads *this from a binary file.")
      //.def(SerializableVisitor<polynomial3_t>())
    ;

    /** END polynomial function**/
    /** BEGIN piecewise curve function **/
    class_<piecewise_polynomial3_curve_t>
    ("piecewise_polynomial3_curve", init<>())
      .def("__init__", make_constructor(&wrapPiecewisePolynomial3CurveConstructor),
           "Create a peicewise-polynomial curve containing the given polynomial curve.")
      .def("min", &piecewise_polynomial3_curve_t::min,"Set the LOWER bound on interval definition of the curve.")
      .def("max", &piecewise_polynomial3_curve_t::max,"Set the HIGHER bound on interval definition of the curve.")
      .def("dim", &piecewise_polynomial3_curve_t::dim)
      .def("__call__", &piecewise_polynomial3_curve_t::operator(),"Evaluate the curve at the given time.")
      .def("derivate", &piecewise_polynomial3_curve_t::derivate,"Evaluate the derivative of order N of curve at time t.",args("self","t","N"))
      .def("add_curve", &piecewise_polynomial3_curve_t::add_curve,
           "Add a new curve to piecewise curve, which should be defined in T_{min},T_{max}] "
           "where T_{min} is equal toT_{max} of the actual piecewise curve.")
      .def("is_continuous", &piecewise_polynomial3_curve_t::is_continuous,"Check if the curve is continuous at the given order.")
      .def("saveAsText", &piecewise_polynomial3_curve_t::saveAsText<piecewise_polynomial3_curve_t>,bp::args("filename"),"Saves *this inside a text file.")
      .def("loadFromText",&piecewise_polynomial3_curve_t::loadFromText<piecewise_polynomial3_curve_t>,bp::args("filename"),"Loads *this from a text file.")
      .def("saveAsXML",&piecewise_polynomial3_curve_t::saveAsXML<piecewise_polynomial3_curve_t>,bp::args("filename","tag_name"),"Saves *this inside a XML file.")
      .def("loadFromXML",&piecewise_polynomial3_curve_t::loadFromXML<piecewise_polynomial3_curve_t>,bp::args("filename","tag_name"),"Loads *this from a XML file.")
      .def("saveAsBinary",&piecewise_polynomial3_curve_t::saveAsBinary<piecewise_polynomial3_curve_t>,bp::args("filename"),"Saves *this inside a binary file.")
      .def("loadFromBinary",&piecewise_polynomial3_curve_t::loadFromBinary<piecewise_polynomial3_curve_t>,bp::args("filename"),"Loads *this from a binary file.")
      //.def(SerializableVisitor<piecewise_polynomial3_curve_t>())
    ;
    class_<piecewise_bezier3_curve_t>
    ("piecewise_bezier3_curve", init<>())
      .def("__init__", make_constructor(&wrapPiecewiseBezier3CurveConstructor))
      .def("min", &piecewise_bezier3_curve_t::min)
      .def("max", &piecewise_bezier3_curve_t::max)
      .def("dim", &piecewise_bezier3_curve_t::dim)
      .def("__call__", &piecewise_bezier3_curve_t::operator())
      .def("derivate", &piecewise_bezier3_curve_t::derivate)
      .def("add_curve", &piecewise_bezier3_curve_t::add_curve)
      .def("is_continuous", &piecewise_bezier3_curve_t::is_continuous)
      .def("saveAsText", &piecewise_bezier3_curve_t::saveAsText<piecewise_bezier3_curve_t>,bp::args("filename"),"Saves *this inside a text file.")
      .def("loadFromText",&piecewise_bezier3_curve_t::loadFromText<piecewise_bezier3_curve_t>,bp::args("filename"),"Loads *this from a text file.")
      .def("saveAsXML",&piecewise_bezier3_curve_t::saveAsXML<piecewise_bezier3_curve_t>,bp::args("filename","tag_name"),"Saves *this inside a XML file.")
      .def("loadFromXML",&piecewise_bezier3_curve_t::loadFromXML<piecewise_bezier3_curve_t>,bp::args("filename","tag_name"),"Loads *this from a XML file.")
      .def("saveAsBinary",&piecewise_bezier3_curve_t::saveAsBinary<piecewise_bezier3_curve_t>,bp::args("filename"),"Saves *this inside a binary file.")
      .def("loadFromBinary",&piecewise_bezier3_curve_t::loadFromBinary<piecewise_bezier3_curve_t>,bp::args("filename"),"Loads *this from a binary file.")
      //.def(SerializableVisitor<piecewise_bezier3_curve_t>())
    ;
    class_<piecewise_cubic_hermite_curve_t>
    ("piecewise_cubic_hermite3_curve", init<>())
      .def("__init__", make_constructor(&wrapPiecewiseCubicHermite3CurveConstructor))
      .def("min", &piecewise_cubic_hermite_curve_t::min)
      .def("max", &piecewise_cubic_hermite_curve_t::max)
      .def("dim", &piecewise_cubic_hermite_curve_t::dim)
      .def("__call__", &piecewise_cubic_hermite_curve_t::operator())
      .def("derivate", &piecewise_cubic_hermite_curve_t::derivate)
      .def("add_curve", &piecewise_cubic_hermite_curve_t::add_curve)
      .def("is_continuous", &piecewise_cubic_hermite_curve_t::is_continuous)
      .def("saveAsText", &piecewise_cubic_hermite_curve_t::saveAsText<piecewise_cubic_hermite_curve_t>,bp::args("filename"),"Saves *this inside a text file.")
      .def("loadFromText",&piecewise_cubic_hermite_curve_t::loadFromText<piecewise_cubic_hermite_curve_t>,bp::args("filename"),"Loads *this from a text file.")
      .def("saveAsXML",&piecewise_cubic_hermite_curve_t::saveAsXML<piecewise_cubic_hermite_curve_t>,bp::args("filename","tag_name"),"Saves *this inside a XML file.")
      .def("loadFromXML",&piecewise_cubic_hermite_curve_t::loadFromXML<piecewise_cubic_hermite_curve_t>,bp::args("filename","tag_name"),"Loads *this from a XML file.")
      .def("saveAsBinary",&piecewise_cubic_hermite_curve_t::saveAsBinary<piecewise_cubic_hermite_curve_t>,bp::args("filename"),"Saves *this inside a binary file.")
      .def("loadFromBinary",&piecewise_cubic_hermite_curve_t::loadFromBinary<piecewise_cubic_hermite_curve_t>,bp::args("filename"),"Loads *this from a binary file.")
      //.def(SerializableVisitor<piecewise_cubic_hermite_curve_t>())
    ;
    /** END piecewise curve function **/
    /** BEGIN exact_cubic curve**/
    class_<exact_cubic3_t>
    ("exact_cubic3", init<>())
      .def("__init__", make_constructor(&wrapExactCubic3Constructor))
      .def("__init__", make_constructor(&wrapExactCubic3ConstructorConstraint))
      .def("min", &exact_cubic3_t::min)
      .def("max", &exact_cubic3_t::max)
      .def("dim", &exact_cubic3_t::dim)
      .def("__call__", &exact_cubic3_t::operator())
      .def("derivate", &exact_cubic3_t::derivate)
      .def("getNumberSplines", &exact_cubic3_t::getNumberSplines)
      .def("getSplineAt", &exact_cubic3_t::getSplineAt)
      .def("saveAsText", &exact_cubic3_t::saveAsText<exact_cubic3_t>,bp::args("filename"),"Saves *this inside a text file.")
      .def("loadFromText",&exact_cubic3_t::loadFromText<exact_cubic3_t>,bp::args("filename"),"Loads *this from a text file.")
      .def("saveAsXML",&exact_cubic3_t::saveAsXML<exact_cubic3_t>,bp::args("filename","tag_name"),"Saves *this inside a XML file.")
      .def("loadFromXML",&exact_cubic3_t::loadFromXML<exact_cubic3_t>,bp::args("filename","tag_name"),"Loads *this from a XML file.")
      .def("saveAsBinary",&exact_cubic3_t::saveAsBinary<exact_cubic3_t>,bp::args("filename"),"Saves *this inside a binary file.")
      .def("loadFromBinary",&exact_cubic3_t::loadFromBinary<exact_cubic3_t>,bp::args("filename"),"Loads *this from a binary file.")
      //.def(SerializableVisitor<exact_cubic3_t>())
    ;
    /** END exact_cubic curve**/
    /** BEGIN cubic_hermite_spline **/
    class_<cubic_hermite_spline_t>
    ("cubic_hermite_spline3", init<>())
      .def("__init__", make_constructor(&wrapCubicHermiteSplineConstructor))
      .def("min", &cubic_hermite_spline_t::min)
      .def("max", &cubic_hermite_spline_t::max)
      .def("dim", &cubic_hermite_spline_t::dim)
      .def("__call__", &cubic_hermite_spline_t::operator())
      .def("derivate", &cubic_hermite_spline_t::derivate)
      .def("saveAsText", &cubic_hermite_spline_t::saveAsText<cubic_hermite_spline_t>,bp::args("filename"),"Saves *this inside a text file.")
      .def("loadFromText",&cubic_hermite_spline_t::loadFromText<cubic_hermite_spline_t>,bp::args("filename"),"Loads *this from a text file.")
      .def("saveAsXML",&cubic_hermite_spline_t::saveAsXML<cubic_hermite_spline_t>,bp::args("filename","tag_name"),"Saves *this inside a XML file.")
      .def("loadFromXML",&cubic_hermite_spline_t::loadFromXML<cubic_hermite_spline_t>,bp::args("filename","tag_name"),"Loads *this from a XML file.")
      .def("saveAsBinary",&cubic_hermite_spline_t::saveAsBinary<cubic_hermite_spline_t>,bp::args("filename"),"Saves *this inside a binary file.")
      .def("loadFromBinary",&cubic_hermite_spline_t::loadFromBinary<cubic_hermite_spline_t>,bp::args("filename"),"Loads *this from a binary file.")
      //.def(SerializableVisitor<cubic_hermite_spline_t>())
    ;
    /** END cubic_hermite_spline **/
    /** BEGIN curve constraints**/
    class_<curve_constraints3_t>
    ("curve_constraints3", init<>())
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
    def("polynomial_from_bezier", polynomial_from_curve<polynomial3_t,bezier3_t>);
    def("polynomial_from_hermite", polynomial_from_curve<polynomial3_t,cubic_hermite_spline_t>);
    def("bezier_from_hermite", bezier_from_curve<bezier3_t,cubic_hermite_spline_t>);
    def("bezier_from_polynomial", bezier_from_curve<bezier3_t,polynomial3_t>);
    def("hermite_from_bezier", hermite_from_curve<cubic_hermite_spline_t, bezier3_t>);
    def("hermite_from_polynomial", hermite_from_curve<cubic_hermite_spline_t, polynomial3_t>);
    /** END curves conversion**/
  } // End BOOST_PYTHON_MODULE
} // namespace curves
