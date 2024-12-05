#ifdef CURVES_WITH_PINOCCHIO_SUPPORT
#include <pinocchio/spatial/motion.hpp>
#include <pinocchio/spatial/se3.hpp>
#endif  // CURVES_WITH_PINOCCHIO_SUPPORT

#include <boost/python.hpp>
#include <boost/python/call_method.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/module_init.hpp>
#include <boost/python/register_ptr_to_python.hpp>
#include <boost/ref.hpp>
#include <boost/utility.hpp>
#include <ndcurves/serialization/curves.hpp>

#include "archive_python_binding.h"
#include "optimization_python.h"
#include "python_variables.h"

namespace ndcurves {
using namespace boost::python;

/* base wrap of curve_abc and others parent abstract class: must implement all
 * pure virtual methods */
struct curve_abc_callback : curve_abc_t {
  curve_abc_callback(PyObject* p) : self(p) {}
  virtual point_t operator()(const real t) const {
    return call_method<point_t>(self, "operator()", t);
  }
  virtual point_t derivate(const real t, const std::size_t n) const {
    return call_method<point_t>(self, "derivate", t, n);
  }
  virtual curve_t* compute_derivate_ptr(const std::size_t n) const {
    return call_method<curve_t*>(self, "compute_derivate", n);
  }
  virtual std::size_t dim() const {
    return call_method<std::size_t>(self, "dim");
  }
  virtual real min() const { return call_method<real>(self, "min"); }
  virtual real max() const { return call_method<real>(self, "max"); }
  virtual std::size_t degree() const {
    return call_method<std::size_t>(self, "degree");
  }
  virtual bool isApprox(
      const curve_t* other,
      const real prec = Eigen::NumTraits<real>::dummy_precision()) const {
    return call_method<bool>(self, "isApprox", other, prec);
  }
  PyObject* self;
};
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(curve_abc_t_isEquivalent_overloads,
                                       curve_abc_t::isEquivalent, 1, 3)

struct curve_3_callback : curve_3_t {
  curve_3_callback(PyObject* p) : self(p) {}
  virtual point3_t operator()(const real t) const {
    return call_method<point3_t>(self, "operator()", t);
  }
  virtual point3_t derivate(const real t, const std::size_t n) const {
    return call_method<point3_t>(self, "derivate", t, n);
  }
  virtual curve_t* compute_derivate_ptr(const std::size_t n) const {
    return call_method<curve_t*>(self, "compute_derivate", n);
  }
  virtual std::size_t dim() const {
    return call_method<std::size_t>(self, "dim");
  }
  virtual real min() const { return call_method<real>(self, "min"); }
  virtual real max() const { return call_method<real>(self, "max"); }
  virtual std::size_t degree() const {
    return call_method<std::size_t>(self, "degree");
  }
  virtual bool isApprox(
      const curve_t* other,
      const real prec = Eigen::NumTraits<real>::dummy_precision()) const {
    return call_method<bool>(self, "isApprox", other, prec);
  }
  PyObject* self;
};
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(curve_3_t_isEquivalent_overloads,
                                       curve_3_t::isEquivalent, 1, 3)

struct curve_rotation_callback : curve_rotation_t {
  curve_rotation_callback(PyObject* p) : self(p) {}
  virtual curve_rotation_t::point_t operator()(const real t) const {
    return call_method<curve_rotation_t::point_t>(self, "operator()", t);
  }
  virtual curve_rotation_t::point_derivate_t derivate(
      const real t, const std::size_t n) const {
    return call_method<curve_rotation_t::point_derivate_t>(self, "derivate", t,
                                                           n);
  }
  virtual curve_rotation_t::curve_derivate_t* compute_derivate_ptr(
      const std::size_t n) const {
    return call_method<curve_rotation_t::curve_derivate_t*>(
        self, "compute_derivate", n);
  }
  virtual std::size_t dim() const {
    return call_method<std::size_t>(self, "dim");
  }
  virtual real min() const { return call_method<real>(self, "min"); }
  virtual real max() const { return call_method<real>(self, "max"); }
  virtual std::size_t degree() const {
    return call_method<std::size_t>(self, "degree");
  }
  virtual bool isApprox(
      const curve_t* other,
      const real prec = Eigen::NumTraits<real>::dummy_precision()) const {
    return call_method<bool>(self, "isApprox", other, prec);
  }
  PyObject* self;
};
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(curve_rotation_t_isEquivalent_overloads,
                                       curve_rotation_t::isEquivalent, 1, 3)

struct curve_SE3_callback : curve_SE3_t {
  curve_SE3_callback(PyObject* p) : self(p) {}
  virtual curve_SE3_t::point_t operator()(const real t) const {
    return call_method<curve_SE3_t::point_t>(self, "operator()", t);
  }
  virtual curve_SE3_t::point_derivate_t derivate(const real t,
                                                 const std::size_t n) const {
    return call_method<curve_SE3_t::point_derivate_t>(self, "derivate", t, n);
  }
  virtual curve_SE3_t::curve_derivate_t* compute_derivate_ptr(
      const std::size_t n) const {
    return call_method<curve_SE3_t::curve_derivate_t*>(self, "compute_derivate",
                                                       n);
  }
  virtual std::size_t dim() const {
    return call_method<std::size_t>(self, "dim");
  }
  virtual real min() const { return call_method<real>(self, "min"); }
  virtual real max() const { return call_method<real>(self, "max"); }
  virtual std::size_t degree() const {
    return call_method<std::size_t>(self, "degree");
  }
  virtual bool isApprox(
      const curve_t* other,
      const real prec = Eigen::NumTraits<real>::dummy_precision()) const {
    return call_method<bool>(self, "isApprox", other, prec);
  }
  PyObject* self;
};
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(curve_SE3_t_isEquivalent_overloads,
                                       curve_SE3_t::isEquivalent, 1, 3)

/* end base wrap of curve_abc */

/* Structure used to define pickle serialization of python curves */
template <typename Curve>
struct curve_pickle_suite : pickle_suite {
  static object getstate(const Curve& curve) {
    std::ostringstream os;
    boost::archive::text_oarchive oa(os);
    oa << curve;
    return str(os.str());
  }

  static void setstate(Curve& curve, object entries) {
    str s = extract<str>(entries)();
    std::string st = extract<std::string>(s)();
    std::istringstream is(st);
    boost::archive::text_iarchive ia(is);
    ia >> curve;
  }
};

template <class C>
struct CopyableVisitor : public def_visitor<CopyableVisitor<C>> {
  template <class PyClass>
  void visit(PyClass& cl) const {
    cl.def("copy", &copy, bp::arg("self"), "Returns a copy of *this.");
    cl.def("__copy__", &copy, bp::arg("self"), "Returns a copy of *this.");
    cl.def("__deepcopy__", &deepcopy, bp::args("self", "memo"),
           "Returns a deep copy of *this.");
  }

 private:
  static C copy(const C& self) { return C(self); }
  static C deepcopy(const C& self, bp::dict) { return C(self); }
};

/* Template constructor bezier */
template <typename Bezier, typename PointList, typename T_Point>
Bezier* wrapBezierConstructorTemplate(const PointList& array,
                                      const real T_min = 0.,
                                      const real T_max = 1.) {
  T_Point asVector = vectorFromEigenArray<PointList, T_Point>(array);
  return new Bezier(asVector.begin(), asVector.end(), T_min, T_max);
}

template <typename Bezier, typename PointList, typename T_Point,
          typename CurveConstraints>
Bezier* wrapBezierConstructorConstraintsTemplate(
    const PointList& array, const CurveConstraints& constraints,
    const real T_min = 0., const real T_max = 1.) {
  T_Point asVector = vectorFromEigenArray<PointList, T_Point>(array);
  return new Bezier(asVector.begin(), asVector.end(), constraints, T_min,
                    T_max);
}
/* End Template constructor bezier */
/* Helper converter constraintsX -> constraints 3 */
curve_constraints3_t convertToConstraints3(curve_constraints_t constraintsX) {
  curve_constraints3_t constraints3(3);
  constraints3.init_vel = point3_t(constraintsX.init_vel);
  constraints3.init_acc = point3_t(constraintsX.init_acc);
  constraints3.end_vel = point3_t(constraintsX.end_vel);
  constraints3.end_acc = point3_t(constraintsX.end_acc);
  return constraints3;
}
/* END helper converter constraintsX -> constraints 3 */

/*3D constructors bezier */
bezier3_t* wrapBezier3Constructor(const pointX_list_t& array) {
  return wrapBezierConstructorTemplate<bezier3_t, pointX_list_t, t_point3_t>(
      array);
}
bezier3_t* wrapBezier3ConstructorBounds(const pointX_list_t& array,
                                        const real T_min, const real T_max) {
  return wrapBezierConstructorTemplate<bezier3_t, pointX_list_t, t_point3_t>(
      array, T_min, T_max);
}
bezier3_t* wrapBezier3ConstructorConstraints(
    const pointX_list_t& array, const curve_constraints_t& constraints) {
  return wrapBezierConstructorConstraintsTemplate<
      bezier3_t, pointX_list_t, t_point3_t, curve_constraints3_t>(
      array, convertToConstraints3(constraints));
}
bezier3_t* wrapBezier3ConstructorBoundsConstraints(
    const pointX_list_t& array, const curve_constraints_t& constraints,
    const real T_min, const real T_max) {
  return wrapBezierConstructorConstraintsTemplate<
      bezier3_t, pointX_list_t, t_point3_t, curve_constraints3_t>(
      array, convertToConstraints3(constraints), T_min, T_max);
}

pointX_list_t wrapBezier3Waypoints(bezier3_t& self) {
  return vectorToEigenArray<bezier3_t::t_point_t, pointX_list_t>(
      self.waypoints());
}
/*END 3D constructors bezier */

/*constructors bezier */
bezier_t* wrapBezierConstructor(const pointX_list_t& array) {
  return wrapBezierConstructorTemplate<bezier_t, pointX_list_t, t_pointX_t>(
      array);
}
bezier_t* wrapBezierConstructorBounds(const pointX_list_t& array,
                                      const real T_min, const real T_max) {
  return wrapBezierConstructorTemplate<bezier_t, pointX_list_t, t_pointX_t>(
      array, T_min, T_max);
}
bezier_t* wrapBezierConstructorConstraints(
    const pointX_list_t& array, const curve_constraints_t& constraints) {
  return wrapBezierConstructorConstraintsTemplate<
      bezier_t, pointX_list_t, t_pointX_t, curve_constraints_t>(array,
                                                                constraints);
}
bezier_t* wrapBezierConstructorBoundsConstraints(
    const pointX_list_t& array, const curve_constraints_t& constraints,
    const real T_min, const real T_max) {
  return wrapBezierConstructorConstraintsTemplate<
      bezier_t, pointX_list_t, t_pointX_t, curve_constraints_t>(
      array, constraints, T_min, T_max);
}
pointX_list_t wrapBezierWaypoints(bezier_t& self) {
  return vectorToEigenArray<bezier_t::t_point_t, pointX_list_t>(
      self.waypoints());
}
/*END constructors bezier */

/* Wrap Cubic hermite spline */
t_pair_pointX_tangent_t getPairsPointTangent(const pointX_list_t& points,
                                             const pointX_list_t& tangents) {
  t_pair_pointX_tangent_t res;
  if (points.size() != tangents.size()) {
    throw std::length_error("size of points and tangents must be the same");
  }
  for (int i = 0; i < points.cols(); ++i) {
    res.push_back(pair_pointX_tangent_t(points.col(i), tangents.col(i)));
  }
  return res;
}

cubic_hermite_spline_t* wrapCubicHermiteSplineConstructor(
    const pointX_list_t& points, const pointX_list_t& tangents,
    const time_waypoints_t& time_pts) {
  t_pair_pointX_tangent_t ppt = getPairsPointTangent(points, tangents);
  std::vector<real> time_control_pts;
  for (int i = 0; i < time_pts.size(); ++i) {
    time_control_pts.push_back(time_pts[i]);
  }
  return new cubic_hermite_spline_t(ppt.begin(), ppt.end(), time_control_pts);
}
/* End wrap Cubic hermite spline */

/* Wrap polynomial */
polynomial_t* wrapPolynomialConstructor1(const coeff_t& array, const real min,
                                         const real max) {
  return new polynomial_t(array, min, max);
}
polynomial_t* wrapPolynomialConstructor2(const coeff_t& array) {
  return new polynomial_t(array, 0., 1.);
}
polynomial_t* wrapPolynomialConstructorFromBoundaryConditionsDegree1(
    const pointX_t& init, const pointX_t& end, const real min, const real max) {
  return new polynomial_t(init, end, min, max);
}
polynomial_t* wrapPolynomialConstructorFromBoundaryConditionsDegree3(
    const pointX_t& init, const pointX_t& d_init, const pointX_t& end,
    const pointX_t& d_end, const real min, const real max) {
  return new polynomial_t(init, d_init, end, d_end, min, max);
}
polynomial_t* wrapPolynomialConstructorFromBoundaryConditionsDegree5(
    const pointX_t& init, const pointX_t& d_init, const pointX_t& dd_init,
    const pointX_t& end, const pointX_t& d_end, const pointX_t& dd_end,
    const real min, const real max) {
  return new polynomial_t(init, d_init, dd_init, end, d_end, dd_end, min, max);
}
static polynomial_t minimumJerk(const pointX_t& init, const pointX_t& end) {
  return polynomial_t::MinimumJerk(init, end);
}
static polynomial_t minimumJerkWithTiming(const pointX_t& init,
                                          const pointX_t& end, const real min,
                                          const real max) {
  return polynomial_t::MinimumJerk(init, end, min, max);
}
/* End wrap polynomial */

/* Wrap piecewise curve */
piecewise_t* wrapPiecewiseCurveConstructor(const curve_ptr_t& curve) {
  return new piecewise_t(curve);
}
piecewise_t* wrapPiecewisePolynomialCurveEmptyConstructor() {
  return new piecewise_t();
}
piecewise_SE3_t* wrapPiecewiseSE3Constructor(const curve_SE3_ptr_t& curve) {
  return new piecewise_SE3_t(curve);
}
piecewise_SE3_t* wrapPiecewiseSE3EmptyConstructor() {
  return new piecewise_SE3_t();
}

typedef bezier_t::piecewise_curve_t piecewise_bezier_t;
piecewise_bezier_t* wrapPiecewiseBezierConstructor(
    const bezier_t::bezier_curve_ptr_t& curve) {
  return new piecewise_bezier_t(curve);
}
piecewise_bezier_t* wrapPiecewiseBezierEmptyConstructor() {
  return new piecewise_bezier_t();
}
typedef bezier_linear_variable_t::piecewise_curve_t piecewise_linear_bezier_t;
piecewise_linear_bezier_t* wrapPiecewiseBezierLinearConstructor(
    const bezier_linear_variable_t::bezier_curve_ptr_t& curve) {
  return new piecewise_linear_bezier_t(curve);
}
piecewise_linear_bezier_t* wrapPiecewiseBezierLinearEmptyConstructor() {
  return new piecewise_linear_bezier_t();
}

static piecewise_t discretPointToPolynomialC0(
    const pointX_list_t& points, const time_waypoints_t& time_points) {
  t_pointX_t points_list =
      vectorFromEigenArray<pointX_list_t, t_pointX_t>(points);
  t_time_t time_points_list =
      vectorFromEigenVector<time_waypoints_t, t_time_t>(time_points);
  return piecewise_t::convert_discrete_points_to_polynomial<polynomial_t>(
      points_list, time_points_list);
}
static piecewise_t discretPointToPolynomialC1(
    const pointX_list_t& points, const pointX_list_t& points_derivative,
    const time_waypoints_t& time_points) {
  t_pointX_t points_list =
      vectorFromEigenArray<pointX_list_t, t_pointX_t>(points);
  t_pointX_t points_derivative_list =
      vectorFromEigenArray<pointX_list_t, t_pointX_t>(points_derivative);
  t_time_t time_points_list =
      vectorFromEigenVector<time_waypoints_t, t_time_t>(time_points);
  return piecewise_t::convert_discrete_points_to_polynomial<polynomial_t>(
      points_list, points_derivative_list, time_points_list);
}
static piecewise_t discretPointToPolynomialC2(
    const pointX_list_t& points, const pointX_list_t& points_derivative,
    const pointX_list_t& points_second_derivative,
    const time_waypoints_t& time_points) {
  t_pointX_t points_list =
      vectorFromEigenArray<pointX_list_t, t_pointX_t>(points);
  t_pointX_t points_derivative_list =
      vectorFromEigenArray<pointX_list_t, t_pointX_t>(points_derivative);
  t_pointX_t points_second_derivative_list =
      vectorFromEigenArray<pointX_list_t, t_pointX_t>(points_second_derivative);

  t_time_t time_points_list =
      vectorFromEigenVector<time_waypoints_t, t_time_t>(time_points);
  return piecewise_t::convert_discrete_points_to_polynomial<polynomial_t>(
      points_list, points_derivative_list, points_second_derivative_list,
      time_points_list);
}

static piecewise_t load_piecewise_from_text_file(const std::string& filename,
                                                 const real dt,
                                                 const std::size_t dim) {
  return piecewise_t::load_piecewise_from_text_file<polynomial_t>(filename, dt,
                                                                  dim);
}

void addFinalPointC0(piecewise_t& self, const pointX_t& end, const real time) {
  if (self.num_curves() == 0)
    throw std::runtime_error(
        "Piecewise append : you need to add at least one curve before using "
        "append(finalPoint) method.");
  if (self.is_continuous(1) && self.num_curves() > 1)
    std::cout << "Warning: by adding this final point to the piecewise curve, "
                 "you loose C1 continuity and only "
                 "guarantee C0 continuity."
              << std::endl;
  curve_ptr_t pol(new polynomial_t(self(self.max()), end, self.max(), time));
  self.add_curve_ptr(pol);
}
void addFinalPointC1(piecewise_t& self, const pointX_t& end,
                     const pointX_t& d_end, const real time) {
  if (self.num_curves() == 0)
    throw std::runtime_error(
        "Piecewise append : you need to add at least one curve before using "
        "append(finalPoint) method.");
  if (self.is_continuous(2) && self.num_curves() > 1)
    std::cout << "Warning: by adding this final point to the piecewise curve, "
                 "you loose C2 continuity and only "
                 "guarantee C1 continuity."
              << std::endl;
  if (!self.is_continuous(1))
    std::cout << "Warning: the current piecewise curve is not C1 continuous."
              << std::endl;
  curve_ptr_t pol(new polynomial_t(self(self.max()),
                                   self.derivate(self.max(), 1), end, d_end,
                                   self.max(), time));
  self.add_curve_ptr(pol);
}
void addFinalPointC2(piecewise_t& self, const pointX_t& end,
                     const pointX_t& d_end, const pointX_t& dd_end,
                     const real time) {
  if (self.num_curves() == 0)
    throw std::runtime_error(
        "Piecewise append : you need to add at least one curve before using "
        "append(finalPoint) method.");
  if (self.is_continuous(3) && self.num_curves() > 1)
    std::cout << "Warning: by adding this final point to the piecewise curve, "
                 "you loose C3 continuity and only "
                 "guarantee C2 continuity."
              << std::endl;
  if (!self.is_continuous(2))
    std::cout << "Warning: the current piecewise curve is not C2 continuous."
              << std::endl;
  curve_ptr_t pol(new polynomial_t(
      self(self.max()), self.derivate(self.max(), 1),
      self.derivate(self.max(), 2), end, d_end, dd_end, self.max(), time));
  self.add_curve_ptr(pol);
}

/* end wrap piecewise polynomial curve */

/* Wrap piecewise3 curve */
piecewise3_t* wrapPiecewise3CurveConstructor(
    const curve_translation_ptr_t& curve) {
  return new piecewise3_t(curve);
}
piecewise3_t* wrapPiecewise3PolynomialCurveEmptyConstructor() {
  return new piecewise3_t();
}
piecewise_SE3_t* wrapPiecewise3SE3Constructor(const curve_SE3_ptr_t& curve) {
  return new piecewise_SE3_t(curve);
}
piecewise_SE3_t* wrapPiecewise3SE3EmptyConstructor() {
  return new piecewise_SE3_t();
}

typedef bezier3_t::piecewise_curve_t piecewise_bezier3_t;
piecewise_bezier3_t* wrapPiecewise3BezierConstructor(
    const bezier3_t::bezier_curve_ptr_t& curve) {
  return new piecewise_bezier3_t(curve);
}
piecewise_bezier3_t* wrapPiecewise3BezierEmptyConstructor() {
  return new piecewise_bezier3_t();
}
typedef bezier_linear_variable_t::piecewise_curve_t piecewise_linear_bezier_t;
piecewise_linear_bezier_t* wrapPiecewise3BezierLinearConstructor(
    const bezier_linear_variable_t::bezier_curve_ptr_t& curve) {
  return new piecewise_linear_bezier_t(curve);
}
piecewise_linear_bezier_t* wrapPiecewise3BezierLinearEmptyConstructor() {
  return new piecewise_linear_bezier_t();
}

static piecewise3_t discretPointsToPolynomial3C0(
    const pointX_list_t& points, const time_waypoints_t& time_points) {
  t_point3_t points_list =
      vectorFromEigenArray<pointX_list_t, t_point3_t>(points);
  t_time_t time_points_list =
      vectorFromEigenVector<time_waypoints_t, t_time_t>(time_points);
  return piecewise3_t::convert_discrete_points_to_polynomial<polynomial3_t>(
      points_list, time_points_list);
}
static piecewise3_t discretPointsToPolynomial3C1(
    const pointX_list_t& points, const pointX_list_t& points_derivative,
    const time_waypoints_t& time_points) {
  t_point3_t points_list =
      vectorFromEigenArray<pointX_list_t, t_point3_t>(points);
  t_point3_t points_derivative_list =
      vectorFromEigenArray<pointX_list_t, t_point3_t>(points_derivative);
  t_time_t time_points_list =
      vectorFromEigenVector<time_waypoints_t, t_time_t>(time_points);
  return piecewise3_t::convert_discrete_points_to_polynomial<polynomial3_t>(
      points_list, points_derivative_list, time_points_list);
}
static piecewise3_t discretPointsToPolynomial3C2(
    const pointX_list_t& points, const pointX_list_t& points_derivative,
    const pointX_list_t& points_second_derivative,
    const time_waypoints_t& time_points) {
  t_point3_t points_list =
      vectorFromEigenArray<pointX_list_t, t_point3_t>(points);
  t_point3_t points_derivative_list =
      vectorFromEigenArray<pointX_list_t, t_point3_t>(points_derivative);
  t_point3_t points_second_derivative_list =
      vectorFromEigenArray<pointX_list_t, t_point3_t>(points_second_derivative);

  t_time_t time_points_list =
      vectorFromEigenVector<time_waypoints_t, t_time_t>(time_points);
  return piecewise3_t::convert_discrete_points_to_polynomial<polynomial3_t>(
      points_list, points_derivative_list, points_second_derivative_list,
      time_points_list);
}

static piecewise3_t load_piecewise3_from_text_file(const std::string& filename,
                                                   const real dt,
                                                   const std::size_t dim) {
  return piecewise3_t::load_piecewise_from_text_file<polynomial3_t>(filename,
                                                                    dt, dim);
}

void addFinalPoint3C0(piecewise3_t& self, const pointX_t& end,
                      const real time) {
  if (self.num_curves() == 0)
    throw std::runtime_error(
        "Piecewise append : you need to add at least one curve before using "
        "append(finalPoint) method.");
  if (self.is_continuous(1) && self.num_curves() > 1)
    std::cout << "Warning: by adding this final point to the piecewise curve, "
                 "you loose C1 continuity and only "
                 "guarantee C0 continuity."
              << std::endl;
  curve_translation_ptr_t pol(
      new polynomial3_t(self(self.max()), end, self.max(), time));
  self.add_curve_ptr(pol);
}
void addFinalPoint3C1(piecewise3_t& self, const pointX_t& end,
                      const pointX_t& d_end, const real time) {
  if (self.num_curves() == 0)
    throw std::runtime_error(
        "Piecewise append : you need to add at least one curve before using "
        "append(finalPoint) method.");
  if (self.is_continuous(2) && self.num_curves() > 1)
    std::cout << "Warning: by adding this final point to the piecewise curve, "
                 "you loose C2 continuity and only "
                 "guarantee C1 continuity."
              << std::endl;
  if (!self.is_continuous(1))
    std::cout << "Warning: the current piecewise curve is not C1 continuous."
              << std::endl;
  curve_translation_ptr_t pol(new polynomial3_t(self(self.max()),
                                                self.derivate(self.max(), 1),
                                                end, d_end, self.max(), time));
  self.add_curve_ptr(pol);
}
void addFinalPoint3C2(piecewise3_t& self, const pointX_t& end,
                      const pointX_t& d_end, const pointX_t& dd_end,
                      const real time) {
  if (self.num_curves() == 0)
    throw std::runtime_error(
        "Piecewise append : you need to add at least one curve before using "
        "append(finalPoint) method.");
  if (self.is_continuous(3) && self.num_curves() > 1)
    std::cout << "Warning: by adding this final point to the piecewise curve, "
                 "you loose C3 continuity and only "
                 "guarantee C2 continuity."
              << std::endl;
  if (!self.is_continuous(2))
    std::cout << "Warning: the current piecewise curve is not C2 continuous."
              << std::endl;
  curve_translation_ptr_t pol(new polynomial3_t(
      self(self.max()), self.derivate(self.max(), 1),
      self.derivate(self.max(), 2), end, d_end, dd_end, self.max(), time));
  self.add_curve_ptr(pol);
}

/* end wrap piecewise polynomial3 curve */

/* Wrap exact cubic spline */
t_waypoint_t getWayPoints(const coeff_t& array,
                          const time_waypoints_t& time_wp) {
  t_waypoint_t res;
  for (int i = 0; i < array.cols(); ++i) {
    res.push_back(std::make_pair(time_wp(i), array.col(i)));
  }
  return res;
}

exact_cubic_t* wrapExactCubicConstructor(const coeff_t& array,
                                         const time_waypoints_t& time_wp) {
  t_waypoint_t wps = getWayPoints(array, time_wp);
  return new exact_cubic_t(wps.begin(), wps.end());
}

exact_cubic_t* wrapExactCubicConstructorConstraint(
    const coeff_t& array, const time_waypoints_t& time_wp,
    const curve_constraints_t& constraints) {
  t_waypoint_t wps = getWayPoints(array, time_wp);
  return new exact_cubic_t(wps.begin(), wps.end(), constraints);
}

/// For constraints XD
pointX_t get_init_vel(const curve_constraints_t& c) { return c.init_vel; }

pointX_t get_init_acc(const curve_constraints_t& c) { return c.init_acc; }

pointX_t get_init_jerk(const curve_constraints_t& c) { return c.init_jerk; }

pointX_t get_end_vel(const curve_constraints_t& c) { return c.end_vel; }

pointX_t get_end_acc(const curve_constraints_t& c) { return c.end_acc; }

pointX_t get_end_jerk(const curve_constraints_t& c) { return c.end_jerk; }

void set_init_vel(curve_constraints_t& c, const pointX_t& val) {
  c.init_vel = val;
}

void set_init_acc(curve_constraints_t& c, const pointX_t& val) {
  c.init_acc = val;
}

void set_init_jerk(curve_constraints_t& c, const pointX_t& val) {
  c.init_jerk = val;
}

void set_end_vel(curve_constraints_t& c, const pointX_t& val) {
  c.end_vel = val;
}

void set_end_acc(curve_constraints_t& c, const pointX_t& val) {
  c.end_acc = val;
}

void set_end_jerk(curve_constraints_t& c, const pointX_t& val) {
  c.end_jerk = val;
}

bezier_t* bezier_t_compute_primitive_init(const bezier_t* b,
                                          const std::size_t order,
                                          const bezier_t::point_t& init) {
  return new bezier_t(b->compute_primitive(order, init));
}
bezier_t* bezier_t_compute_primitive_zero(const bezier_t* b,
                                          const std::size_t order) {
  return new bezier_t(b->compute_primitive(order));
}

bezier3_t* bezier3_t_compute_primitive_init(const bezier3_t* b,
                                            const std::size_t order,
                                            const bezier3_t::point_t& init) {
  return new bezier3_t(b->compute_primitive(order, init));
}
bezier3_t* bezier3_t_compute_primitive_zero(const bezier3_t* b,
                                            const std::size_t order) {
  return new bezier3_t(b->compute_primitive(order));
}

bezier_linear_variable_t* bezier_linear_variable_t_compute_primitive_init(
    const bezier_linear_variable_t* b, const std::size_t order,
    const bezier_linear_variable_t::point_t* init) {
  return new bezier_linear_variable_t(b->compute_primitive(order, *init));
}
bezier_linear_variable_t* bezier_linear_variable_t_compute_primitive_zero(
    const bezier_linear_variable_t* b, const std::size_t order) {
  return new bezier_linear_variable_t(b->compute_primitive(order));
}

bezier_t* bezier_linear_variable_t_evaluate(const bezier_linear_variable_t* b,
                                            const pointX_t& x) {
  return new bezier_t(
      evaluateLinear<bezier_t, bezier_linear_variable_t>(*b, x));
}

bezier_t::piecewise_curve_t (bezier_t::*splitspe)(
    const bezier_t::vector_x_t&) const = &bezier_t::split;
bezier_linear_variable_t::piecewise_curve_t (
    bezier_linear_variable_t::*split_py)(
    const bezier_linear_variable_t::vector_x_t&) const =
    &bezier_linear_variable_t::split;

/* End wrap exact cubic spline */

/* Wrap SO3Linear */
SO3Linear_t* wrapSO3LinearConstructorFromQuaternion(
    const quaternion_t& init_rot, const quaternion_t& end_rot, const real min,
    const real max) {
  return new SO3Linear_t(init_rot, end_rot, min, max);
}

SO3Linear_t* wrapSO3LinearConstructorFromMatrix(const matrix3_t& init_rot,
                                                const matrix3_t& end_rot,
                                                const real min,
                                                const real max) {
  return new SO3Linear_t(init_rot, end_rot, min, max);
}

/* End wrap SO3Linear */

/* Wrap SE3Curves */

matrix4_t se3Return(const curve_SE3_t& curve, const real t) {
  return curve(t).matrix();
}

matrix3_t se3returnRotation(const curve_SE3_t& curve, const real t) {
  return curve(t).rotation();
}

pointX_t se3returnTranslation(const curve_SE3_t& curve, const real t) {
  return pointX_t(curve(t).translation());
}

SE3Curve_t* wrapSE3CurveFromTransform(const matrix4_t& init_pose,
                                      const matrix4_t& end_pose, const real min,
                                      const real max) {
  return new SE3Curve_t(transform_t(init_pose), transform_t(end_pose), min,
                        max);
}

SE3Curve_t* wrapSE3CurveFromPosAndRotation(const pointX_t& init_pos,
                                           const pointX_t& end_pos,
                                           const matrix3_t& init_rot,
                                           const matrix3_t& end_rot,
                                           const real& t_min,
                                           const real& t_max) {
  return new SE3Curve_t(init_pos, end_pos, init_rot, end_rot, t_min, t_max);
}

SE3Curve_t* wrapSE3CurveFromBezier3Translation(bezier3_t& translation_curve,
                                               const matrix3_t& init_rot,
                                               const matrix3_t& end_rot) {
  std::shared_ptr<bezier3_t> translation = std::make_shared<bezier3_t>(
      translation_curve.waypoints().begin(),
      translation_curve.waypoints().end(), translation_curve.min(),
      translation_curve.max());
  return new SE3Curve_t(translation, init_rot, end_rot);
}

SE3Curve_t* wrapSE3CurveFromTranslation(
    const curve_translation_ptr_t& translation_curve, const matrix3_t& init_rot,
    const matrix3_t& end_rot) {
  return new SE3Curve_t(translation_curve, init_rot, end_rot);
}

SE3Curve_t* wrapSE3CurveFromTwoCurves(
    const curve_translation_ptr_t& translation_curve,
    const curve_rotation_ptr_t& rotation_curve) {
  return new SE3Curve_t(translation_curve, rotation_curve);
}

#ifdef CURVES_WITH_PINOCCHIO_SUPPORT
typedef pinocchio::SE3Tpl<real, 0> SE3_t;
typedef pinocchio::MotionTpl<real, 0> Motion_t;

SE3Curve_t* wrapSE3CurveFromSE3Pinocchio(const SE3_t& init_pose,
                                         const SE3_t& end_pose, const real min,
                                         const real max) {
  return new SE3Curve_t(transform_t(init_pose.toHomogeneousMatrix()),
                        transform_t(end_pose.toHomogeneousMatrix()), min, max);
}

SE3_t se3ReturnPinocchio(const curve_SE3_t& curve, const real t) {
  return SE3_t(curve(t).matrix());
}

Motion_t se3ReturnDerivatePinocchio(const curve_SE3_t& curve, const real t,
                                    const std::size_t order) {
  return Motion_t(curve.derivate(t, order));
}
#endif  // CURVES_WITH_PINOCCHIO_SUPPORT
/* End wrap SE3Curves */

/* Wrap piecewiseSE3Curves */
#ifdef CURVES_WITH_PINOCCHIO_SUPPORT
typedef pinocchio::SE3Tpl<real, 0> SE3_t;
typedef pinocchio::MotionTpl<real, 0> Motion_t;

void addFinalSE3(piecewise_SE3_t& self, const SE3_t& end, const real time) {
  if (self.num_curves() == 0)
    throw std::runtime_error(
        "Piecewise append : you need to add at least one curve before using "
        "append(finalPoint) method.");
  if (self.is_continuous(1) && self.num_curves() > 1)
    std::cout << "Warning: by adding this final transform to the piecewise "
                 "curve, you loose C1 continuity and only "
                 "guarantee C0 continuity."
              << std::endl;
  SE3Curve_t curve(self(self.max()), transform_t(end.toHomogeneousMatrix()),
                   self.max(), time);
  self.add_curve(curve);
}

#endif  // CURVES_WITH_PINOCCHIO_SUPPORT

void addFinalTransform(piecewise_SE3_t& self, const matrix4_t& end,
                       const real time) {
  if (self.num_curves() == 0)
    throw std::runtime_error(
        "Piecewise append : you need to add at least one curve before using "
        "append(finalPoint) method.");
  if (self.is_continuous(1) && self.num_curves() > 1)
    std::cout << "Warning: by adding this final transform to the piecewise "
                 "curve, you loose C1 continuity and only "
                 "guarantee C0 continuity."
              << std::endl;
  SE3Curve_t curve(self(self.max()), transform_t(end), self.max(), time);
  self.add_curve(curve);
}

/* End wrap piecewiseSE3Curves */

/* Wrap constant */
constant_t* wrapConstantConstructorTime(const pointX_t& value, const real min,
                                        const real max) {
  return new constant_t(value, min, max);
}
constant_t* wrapConstantConstructor(const pointX_t& value) {
  return new constant_t(value);
}
/* End wrap constant */
/* Wrap constant 3*/
constant3_t* wrapConstant3ConstructorTime(const pointX_t& value, const real min,
                                          const real max) {
  return new constant3_t(value, min, max);
}
constant3_t* wrapConstant3Constructor(const pointX_t& value) {
  return new constant3_t(value);
}
/* End wrap constant 3*/
/* Wrap sinusoidal */
sinusoidal_t* wrapSinusoidalConstructorTime(const pointX_t& p0,
                                            const pointX_t& amplitude,
                                            const real T, const real phi,
                                            const real min, const real max) {
  return new sinusoidal_t(p0, amplitude, T, phi, min, max);
}
sinusoidal_t* wrapSinusoidalConstructor(const pointX_t& p0,
                                        const pointX_t& amplitude, const real T,
                                        const real phi) {
  return new sinusoidal_t(p0, amplitude, T, phi);
}
sinusoidal_t* wrapSinusoidalConstructorStationaryTime(const real time_traj,
                                                      const pointX_t& p_init,
                                                      const pointX_t& p_final,
                                                      const real min,
                                                      const real max) {
  return new sinusoidal_t(time_traj, p_init, p_final, min, max);
}
sinusoidal_t* wrapSinusoidalConstructorStationary(const real time_traj,
                                                  const pointX_t& p_init,
                                                  const pointX_t& p_final) {
  return new sinusoidal_t(time_traj, p_init, p_final);
}
/* End wrap sinusoidal */
// TO DO : Replace all load and save function for serialization in class by
// using
//         SerializableVisitor in archive_python_binding.
BOOST_PYTHON_MODULE(ndcurves) {
  /** BEGIN eigenpy init**/
  eigenpy::enableEigenPy();
  ENABLE_SPECIFIC_MATRIX_TYPE(pointX_t);
  ENABLE_SPECIFIC_MATRIX_TYPE(point3_t);
  ENABLE_SPECIFIC_MATRIX_TYPE(point6_t);
  ENABLE_SPECIFIC_MATRIX_TYPE(pointX_list_t);
  ENABLE_SPECIFIC_MATRIX_TYPE(coeff_t);
  ENABLE_SPECIFIC_MATRIX_TYPE(matrix3_t);
  ENABLE_SPECIFIC_MATRIX_TYPE(matrix4_t);
  // ENABLE_SPECIFIC_MATRIX_TYPE(quaternion_t);
  eigenpy::exposeQuaternion();
  /*eigenpy::exposeAngleAxis();
  eigenpy::exposeQuaternion();*/
  /** END eigenpy init**/
  /** Expose base abstracts class for each dimension/type : **/
  class_<curve_abc_t, boost::noncopyable, std::shared_ptr<curve_abc_callback>>(
      "curve")
      .def("__call__", &curve_abc_t::operator(),
           "Evaluate the curve at the given time.", args("self", "t"))
      .def("derivate", &curve_abc_t::derivate,
           "Evaluate the derivative of order N of curve at time t.",
           args("self", "t", "N"))
      .def("isEquivalent", &curve_abc_t::isEquivalent,
           curve_abc_t_isEquivalent_overloads(
               (bp::arg("other"),
                bp::arg("prec") = Eigen::NumTraits<double>::dummy_precision(),
                bp::arg("order") = 5),
               "isEquivalent check if self and other are approximately equal "
               "by values, given a precision threshold."))
      .def("compute_derivate", &curve_abc_t::compute_derivate_ptr,
           return_value_policy<manage_new_object>(),
           "Return the derivative of *this at the order N.", args("self", "N"))
      .def("min", &curve_abc_t::min,
           "Get the LOWER bound on interval definition of the curve.")
      .def("max", &curve_abc_t::max,
           "Get the HIGHER bound on interval definition of the curve.")
      .def("dim", &curve_abc_t::dim, "Get the dimension of the curve.")
      .def("degree", &curve_abc_t::degree,
           "Get the degree of the representation of the curve (if applicable).")
      .def("saveAsText", pure_virtual(&curve_abc_t::saveAsText<curve_abc_t>),
           bp::args("filename"), "Saves *this inside a text file.")
      .def("loadFromText",
           pure_virtual(&curve_abc_t::loadFromText<curve_abc_t>),
           bp::args("filename"), "Loads *this from a text file.")
      .def("saveAsXML", pure_virtual(&curve_abc_t::saveAsXML<curve_abc_t>),
           bp::args("filename", "tag_name"), "Saves *this inside a XML file.")
      .def("loadFromXML", pure_virtual(&curve_abc_t::loadFromXML<curve_abc_t>),
           bp::args("filename", "tag_name"), "Loads *this from a XML file.")
      .def("saveAsBinary",
           pure_virtual(&curve_abc_t::saveAsBinary<curve_abc_t>),
           bp::args("filename"), "Saves *this inside a binary file.")
      .def("loadFromBinary",
           pure_virtual(&curve_abc_t::loadFromBinary<curve_abc_t>),
           bp::args("filename"), "Loads *this from a binary file.")
      .def_pickle(curve_pickle_suite<curve_abc_t>());

  class_<curve_3_t, boost::noncopyable, bases<curve_abc_t>,
         std::shared_ptr<curve_3_callback>>("curve3")
      .def("__call__", &curve_3_t::operator(),
           "Evaluate the curve at the given time.", args("self", "t"))
      .def("derivate", &curve_3_t::derivate,
           "Evaluate the derivative of order N of curve at time t.",
           args("self", "t", "N"))
      .def("isEquivalent", &curve_3_t::isEquivalent,
           curve_3_t_isEquivalent_overloads(
               (bp::arg("other"),
                bp::arg("prec") = Eigen::NumTraits<double>::dummy_precision(),
                bp::arg("order") = 5),
               "isEquivalent check if self and other are approximately equal "
               "by values, given a precision threshold."))
      .def("compute_derivate", &curve_3_t::compute_derivate_ptr,
           return_value_policy<manage_new_object>(),
           "Return the derivative of *this at the order N.", args("self", "N"))
      .def("min", &curve_3_t::min,
           "Get the LOWER bound on interval definition of the curve.")
      .def("max", &curve_3_t::max,
           "Get the HIGHER bound on interval definition of the curve.")
      .def("dim", &curve_3_t::dim, "Get the dimension of the curve.")
      .def("degree", &curve_3_t::degree,
           "Get the degree of the representation of the curve (if applicable).")
      .def_pickle(curve_pickle_suite<curve_3_t>());

  class_<curve_rotation_t, boost::noncopyable, bases<curve_abc_t>,
         std::shared_ptr<curve_rotation_callback>>("curve_rotation")
      .def("__call__", &curve_rotation_t::operator(),
           "Evaluate the curve at the given time.", args("self", "t"))
      .def("derivate", &curve_rotation_t::derivate,
           "Evaluate the derivative of order N of curve at time t.",
           args("self", "t", "N"))
      .def("isEquivalent", &curve_rotation_t::isEquivalent,
           curve_rotation_t_isEquivalent_overloads(
               (bp::arg("other"),
                bp::arg("prec") = Eigen::NumTraits<double>::dummy_precision(),
                bp::arg("order") = 5),
               "isEquivalent check if self and other are approximately equal "
               "by values, given a precision threshold."))
      .def("compute_derivate", &curve_rotation_t::compute_derivate_ptr,
           return_value_policy<manage_new_object>(),
           "Return the derivative of *this at the order N.", args("self", "N"))
      .def("min", &curve_rotation_t::min,
           "Get the LOWER bound on interval definition of the curve.")
      .def("max", &curve_rotation_t::max,
           "Get the HIGHER bound on interval definition of the curve.")
      .def("dim", &curve_rotation_t::dim, "Get the dimension of the curve.")
      .def("degree", &curve_rotation_t::degree,
           "Get the degree of the representation of the curve (if applicable).")
      .def_pickle(curve_pickle_suite<curve_rotation_t>());

  class_<curve_SE3_t, boost::noncopyable, bases<curve_abc_t>,
         std::shared_ptr<curve_SE3_callback>>("curve_SE3")
      .def("__call__", &se3Return,
           "Evaluate the curve at the given time. Return as an homogeneous "
           "matrix.",
           args("self", "t"))
      .def("derivate", &curve_SE3_t::derivate,
           "Evaluate the derivative of order N of curve at time t. Return as a "
           "vector 6.",
           args("self", "t", "N"))
      .def("isEquivalent", &curve_SE3_t::isEquivalent,
           curve_SE3_t_isEquivalent_overloads(
               (bp::arg("other"),
                bp::arg("prec") = Eigen::NumTraits<double>::dummy_precision(),
                bp::arg("order") = 5),
               "isEquivalent check if self and other are approximately equal "
               "by values, given a precision threshold."))
      .def("compute_derivate", &curve_SE3_t::compute_derivate_ptr,
           return_value_policy<manage_new_object>(),
           "Return the derivative of *this at the order N.", args("self", "N"))
      .def("min", &curve_SE3_t::min,
           "Get the LOWER bound on interval definition of the curve.")
      .def("max", &curve_SE3_t::max,
           "Get the HIGHER bound on interval definition of the curve.")
      .def("dim", &curve_SE3_t::dim, "Get the dimension of the curve.")
      .def("rotation", &se3returnRotation,
           "Output the rotation (as a 3x3 matrix) at the given time.",
           args("self", "time"))
      .def("translation", &se3returnTranslation,
           "Output the rotation (as a vector 3) at the given time.",
           args("self", "time"))
      .def_pickle(curve_pickle_suite<curve_SE3_t>())
#ifdef CURVES_WITH_PINOCCHIO_SUPPORT
      .def("evaluateAsSE3", &se3ReturnPinocchio,
           "Evaluate the curve at the given time. Return as a pinocchio.SE3 "
           "object",
           args("self", "t"))
      .def("derivateAsMotion", &se3ReturnDerivatePinocchio,
           "Evaluate the derivative of order N of curve at time t. Return as a "
           "pinocchio.Motion",
           args("self", "t", "N"))
#endif  // CURVES_WITH_PINOCCHIO_SUPPORT
      ;
  register_ptr_to_python<curve_ptr_t>();
  register_ptr_to_python<curve3_ptr_t>();
  register_ptr_to_python<curve_translation_ptr_t>();
  register_ptr_to_python<curve_rotation_ptr_t>();
  register_ptr_to_python<curve_SE3_ptr_t>();
  /** END base abstracts class for each dimension/type : **/

  /** BEGIN bezier3 curve**/
  bezier3_t (bezier3_t::*cross_bez3)(const bezier3_t&) const =
      &bezier3_t::cross;
  bezier3_t (bezier3_t::*cross_pointBez3)(const bezier3_t::point_t&) const =
      &bezier3_t::cross;
  class_<bezier3_t, bases<curve_3_t>, std::shared_ptr<bezier3_t>>("bezier3",
                                                                  init<>())
      .def("__init__", make_constructor(&wrapBezier3Constructor))
      .def("__init__", make_constructor(&wrapBezier3ConstructorBounds))
      .def("__init__", make_constructor(&wrapBezier3ConstructorConstraints))
      .def("__init__",
           make_constructor(&wrapBezier3ConstructorBoundsConstraints))
      .def("compute_primitive", &bezier3_t_compute_primitive_init,
           return_value_policy<manage_new_object>())
      .def("compute_primitive", &bezier3_t_compute_primitive_zero,
           return_value_policy<manage_new_object>())
      .def("compute_derivate", &bezier3_t::compute_derivate_ptr,
           return_value_policy<manage_new_object>())
      .def("waypointAtIndex", &bezier3_t::waypointAtIndex)
      .def("waypoints", &wrapBezier3Waypoints)
      .def("elevate", &bezier3_t::elevate, bp::args("order"),
           "Computes a Bezier curve of order degrees higher than the current "
           "curve, but strictly equivalent.")
      .def("elevateSelf", &bezier3_t::elevate_self, bp::args("order"),
           "Elevate the Bezier curve of order degrees higher than the current "
           "curve, but strictly equivalent.")
      .def_readonly("degree", &bezier3_t::degree_)
      .def_readonly("nbWaypoints", &bezier3_t::size_)
      .def("saveAsText", &bezier3_t::saveAsText<bezier3_t>,
           bp::args("filename"), "Saves *this inside a text file.")
      .def("loadFromText", &bezier3_t::loadFromText<bezier3_t>,
           bp::args("filename"), "Loads *this from a text file.")
      .def("saveAsXML", &bezier3_t::saveAsXML<bezier3_t>,
           bp::args("filename", "tag_name"), "Saves *this inside a XML file.")
      .def("loadFromXML", &bezier3_t::loadFromXML<bezier3_t>,
           bp::args("filename", "tag_name"), "Loads *this from a XML file.")
      .def("saveAsBinary", &bezier3_t::saveAsBinary<bezier3_t>,
           bp::args("filename"), "Saves *this inside a binary file.")
      .def("loadFromBinary", &bezier3_t::loadFromBinary<bezier3_t>,
           bp::args("filename"), "Loads *this from a binary file.")
      //.def(SerializableVisitor<bezier_t>())
      .def(bp::self == bp::self)
      .def(bp::self != bp::self)
      .def("cross", cross_bez3, bp::args("other"),
           "Compute the cross product of the current bezier by another bezier. "
           "The cross product p1Xp2 of 2 "
           "polynomials p1 and p2 is defined such that forall t, p1Xp2(t) = "
           "p1(t) X p2(t), with X designing the cross "
           "product. This method of course only makes sense for dimension 3 "
           "curves.")
      .def("cross", cross_pointBez3, bp::args("point"),
           "Compute the cross product PXpt of the current Bezier P by a point "
           "pt, such that for all t, PXpt(t) = P(t) "
           "X pt")
      .def(self *= double())
      .def(self /= double())
      .def(self + bezier3_t())
      .def(self - bezier3_t())
      .def(self += bezier3_t())
      .def(self -= bezier3_t())
      .def(self + bezier3_t::point_t())
      .def(self - bezier3_t::point_t())
      .def(self += bezier3_t::point_t())
      .def(self -= bezier3_t::point_t())
      .def(-self)
      .def(self * double())
      .def(self / double())
      .def(CopyableVisitor<bezier3_t>())
      .def_pickle(curve_pickle_suite<bezier3_t>());
  /** END bezier3 curve**/
  /** BEGIN bezier curve**/
  bezier_t (bezier_t::*cross_bez)(const bezier_t&) const = &bezier_t::cross;
  bezier_t (bezier_t::*cross_pointBez)(const bezier_t::point_t&) const =
      &bezier_t::cross;
  class_<bezier_t, bases<curve_abc_t>, std::shared_ptr<bezier_t>>("bezier",
                                                                  init<>())
      .def("__init__", make_constructor(&wrapBezierConstructor))
      .def("__init__", make_constructor(&wrapBezierConstructorBounds))
      .def("__init__", make_constructor(&wrapBezierConstructorConstraints))
      .def("__init__",
           make_constructor(&wrapBezierConstructorBoundsConstraints))
      .def("compute_primitive", &bezier_t_compute_primitive_init,
           return_value_policy<manage_new_object>())
      .def("compute_primitive", &bezier_t_compute_primitive_zero,
           return_value_policy<manage_new_object>())
      .def("compute_derivate", &bezier_t::compute_derivate)
      .def("waypointAtIndex", &bezier_t::waypointAtIndex)
      .def("waypoints", &wrapBezierWaypoints)
      .def("elevate", &bezier_t::elevate, bp::args("order"),
           "Computes a Bezier curve of order degrees higher than the current "
           "curve, but strictly equivalent.")
      .def("elevateSelf", &bezier_t::elevate_self, bp::args("order"),
           "Elevate the Bezier curve of order degrees higher than the current "
           "curve, but strictly equivalent.")
      .def_readonly("degree", &bezier_t::degree_)
      .def_readonly("nbWaypoints", &bezier_t::size_)
      .def("split", splitspe)
      .def(bp::self == bp::self)
      .def(bp::self != bp::self)
      .def("cross", cross_bez, bp::args("other"),
           "Compute the cross product of the current bezier by another bezier. "
           "The cross product p1Xp2 of 2 "
           "polynomials p1 and p2 is defined such that forall t, p1Xp2(t) = "
           "p1(t) X p2(t), with X designing the cross "
           "product. This method of course only makes sense for dimension 3 "
           "curves.")
      .def("cross", cross_pointBez, bp::args("point"),
           "Compute the cross product PXpt of the current Bezier P by a point "
           "pt, such that for all t, PXpt(t) = P(t) "
           "X pt")
      .def(self += bezier_t())
      .def(self -= bezier_t())
      .def(self + bezier_t())
      .def(self - bezier_t())
      .def(self += bezier_t::point_t())
      .def(self -= bezier_t::point_t())
      .def(self + bezier_t::point_t())
      .def(self - bezier_t::point_t())
      .def(self *= double())
      .def(self /= double())
      .def(-self)
      .def(self * double())
      .def(self / double())
      .def(SerializableVisitor<bezier_t>())
      .def(CopyableVisitor<bezier_t>())
      .def_pickle(curve_pickle_suite<bezier_t>());
  /** END bezier curve**/
  /** BEGIN variable points bezier curve**/
  class_<matrix_pair>("matrix_pair", no_init)
      .def_readonly("A", &matrix_pair::A)
      .def_readonly("b", &matrix_pair::b);

  class_<LinearBezierVector>("bezierVarVector", no_init)
      .def_readonly("size", &LinearBezierVector::size)
      .def("at", &LinearBezierVector::at,
           return_value_policy<manage_new_object>());

  class_<linear_variable_t>("linear_variable", init<>())
      .def(init<linear_variable_t::vector_x_t>())
      .def(init<linear_variable_t::matrix_x_t>())
      .def(init<linear_variable_t::matrix_x_t, linear_variable_t::vector_x_t>())
      .def(init<linear_variable_t::matrix_x_t, linear_variable_t::vector_x_t>())
      .def("__call__", &linear_variable_t::operator())
      .def(self += linear_variable_t())
      .def(self -= linear_variable_t())
      .def(self *= double())
      .def(self /= double())
      .def(self + linear_variable_t())
      .def(self - linear_variable_t())
      .def(-self)
      .def(self * double())
      .def(self / double())
      .def(self * linear_variable_t())
      .def("B", &linear_variable_t::B,
           return_value_policy<copy_const_reference>())
      .def("c", &linear_variable_t::c,
           return_value_policy<copy_const_reference>())
      .def("size", &linear_variable_t::size)
      .def("isZero", &linear_variable_t::isZero)
      .def("norm", &linear_variable_t::norm)
      .def("cross", &linear_variable_t::cross, bp::args("other"),
           "Compute the cross product of the current linear_variable and the "
           "other. Only works for dimension 3");

  bezier_linear_variable_t (bezier_linear_variable_t::*cross_bez_var)(
      const bezier_linear_variable_t&) const = &bezier_linear_variable_t::cross;
  bezier_linear_variable_t (bezier_linear_variable_t::*cross_point_var)(
      const bezier_linear_variable_t::point_t&) const =
      &bezier_linear_variable_t::cross;
  class_<bezier_linear_variable_t, bases<curve_abc_t>,
         std::shared_ptr<bezier_linear_variable_t>>("bezier_linear_variable",
                                                    no_init)
      .def("__init__", make_constructor(&wrapBezierLinearConstructor))
      .def("__init__", make_constructor(&wrapBezierLinearConstructorBounds))
      .def("min", &bezier_linear_variable_t::min)
      .def("max", &bezier_linear_variable_t::max)
      .def("__call__", &bezier_linear_variable_t::operator())
      .def("evaluate", &bezier_linear_variable_t_evaluate,
           bp::return_value_policy<bp::manage_new_object>())
      .def("derivate", &bezier_linear_variable_t::derivate)
      .def("compute_derivate", &bezier_linear_variable_t::compute_derivate_ptr,
           return_value_policy<manage_new_object>())
      .def("compute_primitive",
           &bezier_linear_variable_t_compute_primitive_init,
           return_value_policy<manage_new_object>())
      .def("compute_primitive",
           &bezier_linear_variable_t_compute_primitive_zero,
           return_value_policy<manage_new_object>())
      .def("split", split_py)
      .def("waypoints", &wayPointsToLists,
           return_value_policy<manage_new_object>())
      .def("waypointAtIndex", &bezier_linear_variable_t::waypointAtIndex)
      .def_readonly("degree", &bezier_linear_variable_t::degree_)
      .def_readonly("nbWaypoints", &bezier_linear_variable_t::size_)
      .def("cross", cross_bez_var, bp::args("other"),
           "Compute the cross product of the current Bezier by another Bezier. "
           "The cross product p1Xp2 of 2 "
           "polynomials p1 and p2 is defined such that forall t, p1Xp2(t) = "
           "p1(t) X p2(t), with X designing the cross "
           "product. This method of course only makes sense for dimension 3 "
           "polynomials.")
      .def("cross", cross_point_var, bp::args("point"),
           "Compute the cross product PXpt of the current Bezier P by a point "
           "pt, such that for all t, PXpt(t) = P(t) "
           "X pt")
      .def(bp::self == bp::self)
      .def(bp::self != bp::self)
      .def(self += bezier_linear_variable_t())
      .def(self -= bezier_linear_variable_t())
      .def(self + bezier_linear_variable_t())
      .def(self - bezier_linear_variable_t())
      .def(self += linear_variable_t())
      .def(self -= linear_variable_t())
      .def(self + linear_variable_t())
      .def(self - linear_variable_t())
      .def(self *= double())
      .def(self /= double())
      .def(-self)
      .def(self * double())
      .def(self / double())
      .def(CopyableVisitor<bezier_linear_variable_t>())
      .def_pickle(curve_pickle_suite<bezier_linear_variable_t>());

  class_<quadratic_variable_t>("cost", no_init)
      .add_property("A", &cost_t_quad)
      .add_property("b", &cost_t_linear)
      .add_property("c", &cost_t_constant);

  /** END variable points bezier curve**/
  /** BEGIN polynomial curve function**/
  polynomial_t (polynomial_t::*cross_pol)(const polynomial_t&) const =
      &polynomial_t::cross;
  polynomial_t (polynomial_t::*cross_point)(const polynomial_t::point_t&)
      const = &polynomial_t::cross;

  class_<polynomial_t, bases<curve_abc_t>, std::shared_ptr<polynomial_t>>(
      "polynomial", init<>())
      .def(
          "__init__",
          make_constructor(&wrapPolynomialConstructor1, default_call_policies(),
                           args("coeffs", "min", "max")),
          "Create polynomial spline from an Eigen matrix of coefficient "
          "defined for t in [min,max]."
          " The matrix should contain one coefficient per column, from the "
          "zero order coefficient,up to the highest "
          "order."
          " Spline order is given by the number of the columns -1.")
      .def("__init__",
           make_constructor(&wrapPolynomialConstructor2,
                            default_call_policies(), arg("coeffs")),
           "Create polynomial spline from an Eigen matrix of coefficient "
           "defined for t in [0,1]."
           " The matrix should contain one coefficient per column, from the "
           "zero order coefficient,up to the highest "
           "order."
           " Spline order is given by the number of the columns -1.")
      .def("__init__",
           make_constructor(
               &wrapPolynomialConstructorFromBoundaryConditionsDegree1,
               default_call_policies(), args("init", "end", "min", "max")),
           "Create a polynomial of degree 1 defined for t in [min,max], "
           "such that c(min) == init and c(max) == end.")
      .def("__init__",
           make_constructor(
               &wrapPolynomialConstructorFromBoundaryConditionsDegree3,
               default_call_policies(),
               args("init", "d_init", "end", "d_end", "min", "max")),
           "Create a polynomial of degree 3 defined for t in [min,max], "
           "such that c(min) == init and c(max) == end"
           " dc(min) == d_init and dc(max) == d_end")
      .def("__init__",
           make_constructor(
               &wrapPolynomialConstructorFromBoundaryConditionsDegree5,
               default_call_policies(),
               args("init", "d_init", "dd_init", "end", "d_end", "dd_end",
                    "min", "max")),
           "Create a polynomial of degree 5 defined for t in [min,max], "
           "such that c(min) == init and c(max) == end"
           " dc(min) == d_init and dc(max) == d_end"
           " ddc(min) == dd_init and ddc(max) == dd_end")
      .def("MinimumJerk", &minimumJerk, args("init", "end"),
           "Build a polynomial curve connecting init to end minimizing the "
           "time integral of the squared jerk,"
           "with a zero initial and final velocity and acceleration."
           "The curve is defined in [0; 1], of duration 1.")
      .def("MinimumJerk", &minimumJerkWithTiming,
           args("init", "end", "t_min", "t_max"),
           "Build a polynomial curve connecting init to end minimizing the "
           "time integral of the squared jerk,"
           "with a zero initial and final velocity and acceleration."
           "The curve is defined in [t_min; t_max], of duration t_max - t_min.")
      .staticmethod("MinimumJerk")
      .def("coeffAtDegree", &polynomial_t::coeffAtDegree)
      .def("coeff", &polynomial_t::coeff)
      .def(SerializableVisitor<polynomial_t>())
      .def("cross", cross_pol, bp::args("other"),
           "Compute the cross product of the current polynomial by another "
           "polynomial. The cross product p1Xp2 of 2 "
           "polynomials p1 and p2 is defined such that forall t, p1Xp2(t) = "
           "p1(t) X p2(t), with X designing the cross "
           "product. This method of course only makes sense for dimension 3 "
           "polynomials.")
      .def("cross", cross_point, bp::args("point"),
           "Compute the cross product PXpt of the current polynomial P by a "
           "point pt, such that for all t, PXpt(t) = "
           "P(t) X pt")
      .def(bp::self == bp::self)
      .def(bp::self != bp::self)
      .def(self += polynomial_t())
      .def(self -= polynomial_t())
      .def(self + polynomial_t())
      .def(self - polynomial_t())
      .def(self += polynomial_t::point_t())
      .def(self -= polynomial_t::point_t())
      .def(self + polynomial_t::point_t())
      .def(self - polynomial_t::point_t())
      .def(self *= double())
      .def(self /= double())
      .def(-self)
      .def(self * double())
      .def(self / double())
      .def(CopyableVisitor<polynomial_t>())
      .def_pickle(curve_pickle_suite<polynomial_t>());

  /** END polynomial function**/
  /** BEGIN piecewise curve function **/
  class_<piecewise_t, bases<curve_abc_t>, std::shared_ptr<piecewise_t>>(
      "piecewise", init<>())
      .def("__init__",
           make_constructor(&wrapPiecewiseCurveConstructor,
                            default_call_policies(), arg("curve")),
           "Create a piecewise curve containing the given curve.")
      .def("FromPointsList", &discretPointToPolynomialC0,
           "Create a piecewise-polynomial connecting exactly all the given "
           "points at the given time. The created "
           "piecewise is C0 continuous.",
           args("points", "time_points"))
      .def("FromPointsList", &discretPointToPolynomialC1,
           "Create a piecewise-polynomial connecting exactly all the given "
           "points at the given time and respect the "
           "given points derivative values. The created piecewise is C1 "
           "continuous.",
           args("points", "points_derivative", "time_points"))
      .def("FromPointsList", &discretPointToPolynomialC2,
           "Create a piecewise-polynomial connecting exactly all the given "
           "points at the given time and respect the "
           "given points derivative and second derivative values. The created "
           "piecewise is C2 continuous.",
           args("points", "points_derivative", "points_second_derivative",
                "time_points"))
      .staticmethod("FromPointsList")
      .def("FromPointsFile", &load_piecewise_from_text_file,
           args("filename", "dt", "dimension"),
           "Create a piecewise-polynomial connecting exactly all the points in "
           "the given text file."
           "The file should contains one points per line, optionally with it's "
           "derivative and second derivatives."
           "Each lines should thus contains dim, 2*dim or 3*dim values")
      .staticmethod("FromPointsFile")
      .def("append", &addFinalPointC0,
           "Append a new polynomial curve of degree 1 at the end of the "
           "piecewise curve, defined between self.max() "
           "and time and connecting exactly self(self.max()) and end",
           args("self", "end", "time"))
      .def("append", &addFinalPointC1,
           "Append a new polynomial curve of degree 3 at the end of the "
           "piecewise curve, defined between self.max() "
           "and time and connecting exactly self(self.max()) and end. It "
           "guarantee C1 continuity and guarantee that "
           "self.derivate(time,1) == d_end",
           args("self", "end", "d_end", "time"))
      .def("append", &addFinalPointC2,
           "Append a new polynomial curve of degree 5 at the end of the "
           "piecewise curve, defined between self.max() "
           "and time and connecting exactly self(self.max()) and end. It "
           "guarantee C2 continuity and guarantee that "
           "self.derivate(time,1) == d_end and self.derivate(time,2) == dd_end",
           args("self", "end", "d_end", "d_end", "time"))
      .def("append", &piecewise_t::add_curve_ptr,
           "Add a new curve to piecewise curve, which should be defined in "
           "T_{min},T_{max}] "
           "where T_{min} is equal toT_{max} of the actual piecewise curve.")
      .def("is_continuous", &piecewise_t::is_continuous,
           "Check if the curve is continuous at the given order.",
           args("self", "order"))
      .def("convert_piecewise_curve_to_polynomial",
           &piecewise_t::convert_piecewise_curve_to_polynomial<polynomial_t>,
           "Convert a piecewise curve to to a piecewise polynomial curve")
      .def("convert_piecewise_curve_to_bezier",
           &piecewise_t::convert_piecewise_curve_to_bezier<bezier_t>,
           "Convert a piecewise curve to to a piecewise bezier curve")
      .def("convert_piecewise_curve_to_cubic_hermite",
           &piecewise_t::convert_piecewise_curve_to_cubic_hermite<
               cubic_hermite_spline_t>,
           "Convert a piecewise curve to to a piecewise cubic hermite spline")
      .def("curve_at_index", &piecewise_t::curve_at_index)
      .def("curve_at_time", &piecewise_t::curve_at_time)
      .def("num_curves", &piecewise_t::num_curves)
      .def(SerializableVisitor<piecewise_t>())
      .def(bp::self == bp::self)
      .def(bp::self != bp::self)
      .def(CopyableVisitor<piecewise_t>())
      .def_pickle(curve_pickle_suite<piecewise_t>());

  class_<piecewise3_t, bases<curve_3_t>, std::shared_ptr<piecewise3_t>>(
      "piecewise3", init<>())
      .def("__init__",
           make_constructor(&wrapPiecewise3CurveConstructor,
                            default_call_policies(), arg("curve")),
           "Create a piecewise curve containing the given curve.")
      .def("FromPointsList", &discretPointsToPolynomial3C0,
           "Create a piecewise-polynomial connecting exactly all the given "
           "points at the given time. The created "
           "piecewise is C0 continuous.",
           args("points", "time_points"))
      .def("FromPointsList", &discretPointsToPolynomial3C1,
           "Create a piecewise-polynomial connecting exactly all the given "
           "points at the given time and respect the "
           "given points derivative values. The created piecewise is C1 "
           "continuous.",
           args("points", "points_derivative", "time_points"))
      .def("FromPointsList", &discretPointsToPolynomial3C2,
           "Create a piecewise-polynomial connecting exactly all the given "
           "points at the given time and respect the "
           "given points derivative and second derivative values. The created "
           "piecewise is C2 continuous.",
           args("points", "points_derivative", "points_second_derivative",
                "time_points"))
      .staticmethod("FromPointsList")
      .def("FromPointsFile", &load_piecewise3_from_text_file,
           args("filename", "dt", "dimension"),
           "Create a piecewise-polynomial connecting exactly all the points in "
           "the given text file."
           "The file should contains one points per line, optionally with it's "
           "derivative and second derivatives."
           "Each lines should thus contains dim, 2*dim or 3*dim values")
      .staticmethod("FromPointsFile")
      .def("append", &addFinalPoint3C0,
           "Append a new polynomial curve of degree 1 at the end of the "
           "piecewise curve, defined between self.max() "
           "and time and connecting exactly self(self.max()) and end",
           args("self", "end", "time"))
      .def("append", &addFinalPoint3C1,
           "Append a new polynomial curve of degree 3 at the end of the "
           "piecewise curve, defined between self.max() "
           "and time and connecting exactly self(self.max()) and end. It "
           "guarantee C1 continuity and guarantee that "
           "self.derivate(time,1) == d_end",
           args("self", "end", "d_end", "time"))
      .def("append", &addFinalPoint3C2,
           "Append a new polynomial curve of degree 5 at the end of the "
           "piecewise curve, defined between self.max() "
           "and time and connecting exactly self(self.max()) and end. It "
           "guarantee C2 continuity and guarantee that "
           "self.derivate(time,1) == d_end and self.derivate(time,2) == dd_end",
           args("self", "end", "d_end", "d_end", "time"))
      .def("append", &piecewise3_t::add_curve_ptr,
           "Add a new curve to piecewise curve, which should be defined in "
           "T_{min},T_{max}] "
           "where T_{min} is equal toT_{max} of the actual piecewise curve.")
      .def("is_continuous", &piecewise3_t::is_continuous,
           "Check if the curve is continuous at the given order.",
           args("self", "order"))
      .def("convert_piecewise_curve_to_polynomial",
           &piecewise3_t::convert_piecewise_curve_to_polynomial<polynomial3_t>,
           "Convert a piecewise curve to to a piecewise polynomial curve")
      .def("convert_piecewise_curve_to_bezier",
           &piecewise3_t::convert_piecewise_curve_to_bezier<bezier3_t>,
           "Convert a piecewise curve to to a piecewise bezier curve")
      .def("curve_at_index", &piecewise3_t::curve_at_index)
      .def("curve_at_time", &piecewise3_t::curve_at_time)
      .def("num_curves", &piecewise3_t::num_curves)
      .def(SerializableVisitor<piecewise3_t>())
      .def(bp::self == bp::self)
      .def(bp::self != bp::self)
      .def(CopyableVisitor<piecewise3_t>())
      .def_pickle(curve_pickle_suite<piecewise3_t>());

  class_<piecewise_bezier_t, bases<curve_abc_t>,
         std::shared_ptr<piecewise_bezier_t>>("piecewise_bezier", init<>())
      .def("__init__",
           make_constructor(&wrapPiecewiseBezierConstructor,
                            default_call_policies(), arg("curve")),
           "Create a peicewise Bezier curve containing the given curve.")
      .def("__init__", make_constructor(&wrapPiecewiseBezierEmptyConstructor),
           "Create an empty piecewise-Beziercurve.")
      .def("append", &piecewise_bezier_t::add_curve_ptr,
           "Add a new curve to piecewise curve, which should be defined in "
           "T_{min},T_{max}] "
           "where T_{min} is equal toT_{max} of the actual piecewise curve.")
      .def("is_continuous", &piecewise_bezier_t::is_continuous,
           "Check if the curve is continuous at the given order.",
           args("self", "order"))
      .def("curve_at_index", &piecewise_bezier_t::curve_at_index)
      .def("curve_at_time", &piecewise_bezier_t::curve_at_time)
      .def("num_curves", &piecewise_bezier_t::num_curves)
      .def(SerializableVisitor<piecewise_bezier_t>())
      .def(bp::self == bp::self)
      .def(bp::self != bp::self)
      .def(CopyableVisitor<piecewise_bezier_t>())
      .def_pickle(curve_pickle_suite<piecewise_bezier_t>());

  class_<piecewise_linear_bezier_t, bases<curve_abc_t>,
         std::shared_ptr<piecewise_linear_bezier_t>>("piecewise_bezier_linear",
                                                     init<>(args("self")))
      .def("__init__",
           make_constructor(&wrapPiecewiseBezierLinearConstructor,
                            default_call_policies(), arg("curve")),
           "Create a peicewise Bezier curve containing the given curve.")
      .def("__init__",
           make_constructor(&wrapPiecewiseBezierLinearEmptyConstructor),
           "Create an empty piecewise-Beziercurve.")
      .def("append", &piecewise_linear_bezier_t::add_curve_ptr,
           "Add a new curve to piecewise curve, which should be defined in "
           "T_{min},T_{max}] "
           "where T_{min} is equal toT_{max} of the actual piecewise curve.")
      .def("is_continuous", &piecewise_linear_bezier_t::is_continuous,
           "Check if the curve is continuous at the given order.",
           args("self", "order"))
      .def("curve_at_index", &piecewise_linear_bezier_t::curve_at_index)
      .def("curve_at_time", &piecewise_linear_bezier_t::curve_at_time)
      .def("num_curves", &piecewise_linear_bezier_t::num_curves)
      .def(SerializableVisitor<piecewise_linear_bezier_t>())
      .def(bp::self == bp::self)
      .def(bp::self != bp::self)
      .def(CopyableVisitor<piecewise_linear_bezier_t>())
      .def_pickle(curve_pickle_suite<piecewise_linear_bezier_t>());

  class_<piecewise_SE3_t, bases<curve_SE3_t>, std::shared_ptr<piecewise_SE3_t>>(
      "piecewise_SE3", init<>())
      .def("__init__",
           make_constructor(&wrapPiecewiseSE3Constructor,
                            default_call_policies(), arg("curve")),
           "Create a piecewise-se3 curve containing the given se3 curve.")
      .def("__init__", make_constructor(&wrapPiecewiseSE3EmptyConstructor),
           "Create an empty piecewise-se3 curve.")
      .def("append", &piecewise_SE3_t::add_curve_ptr,
           "Add a new curve to piecewise curve, which should be defined in "
           "T_{min},T_{max}] "
           "where T_{min} is equal toT_{max} of the actual piecewise curve.",
           args("self", "curve"))
      .def("is_continuous", &piecewise_SE3_t::is_continuous,
           "Check if the curve is continuous at the given order.",
           args("self", "order"))
      .def("curve_at_index", &piecewise_SE3_t::curve_at_index)
      .def("curve_at_time", &piecewise_SE3_t::curve_at_time)
      .def("num_curves", &piecewise_SE3_t::num_curves)
      .def("append", &addFinalTransform,
           "Append a new linear SE3 curve at the end of the piecewise curve, "
           "defined between self.max() "
           "and time and connecting exactly self(self.max()) and end",
           args("self", "end", "time"))
      .def(SerializableVisitor<piecewise_SE3_t>())
      .def(bp::self == bp::self)
      .def(bp::self != bp::self)
      .def(CopyableVisitor<piecewise_SE3_t>())
      .def_pickle(curve_pickle_suite<piecewise_SE3_t>())
#ifdef CURVES_WITH_PINOCCHIO_SUPPORT
      .def("append", &addFinalSE3,
           "Append a new linear SE3 curve at the end of the piecewise curve, "
           "defined between self.max() "
           "and time and connecting exactly self(self.max()) and end",
           args("self", "end", "time"))
#endif  // CURVES_WITH_PINOCCHIO_SUPPORT
      ;

  /** END piecewise curve function **/
  /** BEGIN exact_cubic curve**/
  class_<exact_cubic_t, bases<curve_abc_t>, std::shared_ptr<exact_cubic_t>>(
      "exact_cubic", init<>(args("self")))
      .def("__init__",
           make_constructor(&wrapExactCubicConstructor, default_call_policies(),
                            args("array", "time_wp")))
      .def("__init__",
           make_constructor(&wrapExactCubicConstructorConstraint,
                            default_call_policies(),
                            args("array", "time_wp", "constraints")))
      .def("getNumberSplines", &exact_cubic_t::getNumberSplines, args("self"))
      .def("getSplineAt", &exact_cubic_t::getSplineAt, args("self", "index"))
      .def(SerializableVisitor<exact_cubic_t>())
      .def(bp::self == bp::self)
      .def(bp::self != bp::self)
      .def(CopyableVisitor<exact_cubic_t>())
      .def_pickle(curve_pickle_suite<exact_cubic_t>());

  /** END exact_cubic curve**/
  /** BEGIN cubic_hermite_spline **/
  class_<cubic_hermite_spline_t, bases<curve_abc_t>,
         std::shared_ptr<cubic_hermite_spline_t>>("cubic_hermite_spline",
                                                  init<>(args("self")))
      .def("__init__", make_constructor(&wrapCubicHermiteSplineConstructor,
                                        bp::default_call_policies(),
                                        args("points", "tangents", "times")))
      .def(SerializableVisitor<cubic_hermite_spline_t>())
      .def(bp::self == bp::self)
      .def(bp::self != bp::self)
      .def(CopyableVisitor<cubic_hermite_spline_t>())
      .def(bp::self == bp::self)
      .def_pickle(curve_pickle_suite<cubic_hermite_spline_t>());

  /** END cubic_hermite_spline **/
  /** BEGIN curve constraints**/
  class_<curve_constraints_t>("curve_constraints", init<>())
      .def(bp::init<int>(args("self", "dimension"),
                         "Init with a given dimension."))
      .add_property("init_vel", &get_init_vel, &set_init_vel)
      .add_property("init_acc", &get_init_acc, &set_init_acc)
      .add_property("init_jerk", &get_init_jerk, &set_init_jerk)
      .add_property("end_vel", &get_end_vel, &set_end_vel)
      .add_property("end_acc", &get_end_acc, &set_end_acc)
      .add_property("end_jerk", &get_end_jerk, &set_end_jerk)
      .def("__eq__", &curve_constraints_t::operator==)
      .def("__ne__", &curve_constraints_t::operator!=)
      .def(SerializableVisitor<curve_constraints_t>())
      .def(CopyableVisitor<curve_constraints_t>())
      .def_pickle(curve_pickle_suite<curve_constraints_t>());
  ;
  /** END curve constraints**/
  /** BEGIN bernstein polynomial**/
  class_<bernstein_t>("bernstein",
                      init<const unsigned int, const unsigned int>())
      .def("__call__", &bernstein_t::operator())
      .def(bp::self == bp::self)
      .def(bp::self != bp::self);
  /** END bernstein polynomial**/

  /** BEGIN SO3 Linear**/
  class_<SO3Linear_t, bases<curve_rotation_t>, std::shared_ptr<SO3Linear_t>>(
      "SO3Linear", init<>())
      .def("__init__",
           make_constructor(
               &wrapSO3LinearConstructorFromMatrix, default_call_policies(),
               args("init_rotation", "end_rotation", "min", "max")),
           "Create a SO3 Linear curve between two rotations, defined for t in "
           "[min,max]."
           " The input rotations are expressed as 3x3 matrix.")
      .def("__init__",
           make_constructor(
               &wrapSO3LinearConstructorFromQuaternion, default_call_policies(),
               args("init_rotation", "end_rotation", "min", "max")),
           "Create a SO3 Linear curve between two rotations, defined for t in "
           "[min,max]."
           " The input rotations are expressed as Quaternions.")
      .def("computeAsQuaternion", &SO3Linear_t::computeAsQuaternion,
           "Output the quaternion of the rotation at the given time. This "
           "rotation is obtained by a Spherical Linear "
           "Interpolation between the initial and final rotation.")
      .def(SerializableVisitor<SO3Linear_t>())
      .def(bp::self == bp::self)
      .def(bp::self != bp::self)
      .def(CopyableVisitor<SO3Linear_t>())
      .def_pickle(curve_pickle_suite<SO3Linear_t>());

  /** END  SO3 Linear**/
  /** BEGIN SE3 Curve**/
  class_<SE3Curve_t, bases<curve_SE3_t>, std::shared_ptr<SE3Curve_t>>(
      "SE3Curve", init<>())
      .def("__init__",
           make_constructor(
               &wrapSE3CurveFromTransform, default_call_policies(),
               args("init_transform", "end_transform", "min", "max")),
           "Create a SE3 curve between two transform, defined for t in "
           "[min,max]."
           " Using linear interpolation for translation and slerp for rotation "
           "between init and end."
           " The input transform are expressed as 4x4 matrix.")
      .def("__init__",
           make_constructor(
               &wrapSE3CurveFromPosAndRotation, default_call_policies(),
               args("init_translation", "end_translation", "init_rotation",
                    "end_rotation", "min", "max")),
           "Create a SE3 curve between two transform, defined for t in "
           "[min,max]."
           " Using linear interpolation for translation and slerp for rotation "
           "between init and end."
           " The input translations are expressed as 3d vector and the "
           "rotations as 3x3 matrix.")
      .def("__init__",
           make_constructor(&wrapSE3CurveFromTwoCurves, default_call_policies(),
                            args("translation_curve", "rotation_curve")),
           "Create a SE3 curve from a translation curve and a rotation one."
           "The translation curve should be of dimension 3 and the rotation "
           "one should output 3x3 matrix"
           "Both curves should have the same time bounds.")
      .def("__init__",
           make_constructor(
               &wrapSE3CurveFromTranslation, default_call_policies(),
               args("translation_curve", "init_rotation", "end_rotation")),
           "Create a SE3 curve from a translation curve and two rotation"
           "The translation curve should be of dimension 3, the time "
           "definition of the SE3curve will the same as the "
           "translation curve."
           "The orientation along the SE3Curve will be a slerp between the two "
           "given rotations."
           "The orientations should be represented as 3x3 rotation matrix")
      .def("__init__",
           make_constructor(
               &wrapSE3CurveFromBezier3Translation, default_call_policies(),
               args("translation_curve", "init_rotation", "end_rotation")),
           "Create a SE3 curve from a translation curve and two rotation"
           "The translation curve should be of dimension 3, the time "
           "definition of the SE3curve will the same as the "
           "translation curve."
           "The orientation along the SE3Curve will be a slerp between the two "
           "given rotations."
           "The orientations should be represented as 3x3 rotation matrix")
      .def("translation_curve", &SE3Curve_t::translation_curve,
           "Return a curve corresponding to the translation part of self.")
      .def("rotation_curve", &SE3Curve_t::rotation_curve,
           "Return a curve corresponding to the rotation part of self.")
      .def(SerializableVisitor<SE3Curve_t>())
      .def(bp::self == bp::self)
      .def(bp::self != bp::self)
      .def(CopyableVisitor<SE3Curve_t>())
      .def_pickle(curve_pickle_suite<SE3Curve_t>())
#ifdef CURVES_WITH_PINOCCHIO_SUPPORT
      .def("__init__",
           make_constructor(&wrapSE3CurveFromSE3Pinocchio,
                            default_call_policies(),
                            args("init_SE3", "end_SE3", "min", "max")),
           "Create a SE3 curve between two SE3 objects from Pinocchio, defined "
           "for t in [min,max]."
           " Using linear interpolation for translation and slerp for rotation "
           "between init and end.")
#endif  // CURVES_WITH_PINOCCHIO_SUPPORT
      ;

  /** END SE3 Curve**/
  /** BEGIN constant curve function**/
  class_<constant_t, bases<curve_abc_t>, std::shared_ptr<constant_t>>(
      "constant", init<>())
      .def("__init__",
           make_constructor(&wrapConstantConstructorTime,
                            default_call_policies(),
                            args("value", "min", "max")),
           "Create a constant curve defined for t in [min,max]."
           " This curve always evaluate to the given value and derivate to "
           "zero value.")
      .def("__init__",
           make_constructor(&wrapConstantConstructor, default_call_policies(),
                            arg("value")),
           "Create a constant curve defined for t in [0,inf]."
           " This curve always evaluate to the given value and derivate to "
           "zero value.")
      .def(SerializableVisitor<constant_t>())
      .def(bp::self == bp::self)
      .def(bp::self != bp::self)
      .def(CopyableVisitor<constant_t>())
      .def_pickle(curve_pickle_suite<constant_t>());
  /** END constant function**/
  /** BEGIN constant 3 curve function**/
  class_<constant3_t, bases<curve_3_t>, std::shared_ptr<constant3_t>>(
      "constant3", init<>())
      .def("__init__",
           make_constructor(&wrapConstant3ConstructorTime,
                            default_call_policies(),
                            args("value", "min", "max")),
           "Create a constant curve defined for t in [min,max]."
           " This curve always evaluate to the given value and derivate to "
           "zero value.")
      .def("__init__",
           make_constructor(&wrapConstant3Constructor, default_call_policies(),
                            arg("value")),
           "Create a constant curve defined for t in [0,inf]."
           " This curve always evaluate to the given value and derivate to "
           "zero value.")
      .def(SerializableVisitor<constant3_t>())
      .def(bp::self == bp::self)
      .def(bp::self != bp::self)
      .def(CopyableVisitor<constant3_t>())
      .def_pickle(curve_pickle_suite<constant3_t>());
  /** END constant 3 function**/
  /** BEGIN sinusoidal curve function**/
  class_<sinusoidal_t, bases<curve_abc_t>, std::shared_ptr<sinusoidal_t>>(
      "sinusoidal", init<>())
      .def("__init__",
           make_constructor(&wrapSinusoidalConstructor, default_call_policies(),
                            args("Offset", "Amplitude", "Period", "Phase")),
           "Create a sinusoidal curve defined for t in [0, inf]."
           " c(t) = offset + amplitude * sin(2pi/T * t + phi)")
      .def("__init__",
           make_constructor(
               &wrapSinusoidalConstructorTime, default_call_policies(),
               args("Offset", "Amplitude", "Period", "Phase", "min", "max")),
           "Create a sinusoidal curve defined for t in [min, max]."
           " c(t) = offset + amplitude * sin(2pi/T * t + phi)")
      .def("__init__",
           make_constructor(&wrapSinusoidalConstructorStationary,
                            default_call_policies(),
                            args("duration", "p_init", "p_final")),
           "Create a sinusoidal curve defined for t in [0, inf]."
           "That connect the two stationnary points p_init and p_final in "
           "duration (an half period)")
      .def(
          "__init__",
          make_constructor(&wrapSinusoidalConstructorStationaryTime,
                           default_call_policies(),
                           args("duration", "p_init", "p_final", "min", "max")),
          "Create a sinusoidal curve defined for t in [min, max]."
          "That connect the two stationnary points p_init and p_final in "
          "duration (an half period)")
      .def(SerializableVisitor<sinusoidal_t>())
      .def(bp::self == bp::self)
      .def(bp::self != bp::self)
      .def(CopyableVisitor<sinusoidal_t>())
      .def_pickle(curve_pickle_suite<sinusoidal_t>());
  /** END sinusoidal function**/
  /** BEGIN curves conversion**/
  def("convert_to_polynomial", polynomial_from_curve<polynomial_t>);
  def("convert_to_bezier", bezier_from_curve<bezier_t>);
  def("convert_to_hermite", hermite_from_curve<cubic_hermite_spline_t>);
  /** END curves conversion**/

  optimization::python::exposeOptimization();

#ifdef CURVES_WITH_PINOCCHIO_SUPPORT
  scope().attr("CURVES_WITH_PINOCCHIO_SUPPORT") = true;
#else
  scope().attr("CURVES_WITH_PINOCCHIO_SUPPORT") = false;
#endif

}  // End BOOST_PYTHON_MODULE
}  // namespace ndcurves
