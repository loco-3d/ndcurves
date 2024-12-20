#include <cmath>
#include <ctime>
#include <iostream>
#include <memory>
#include <string>

#include "load_problem.h"
#include "ndcurves/bezier_curve.h"
#include "ndcurves/cubic_hermite_spline.h"
#include "ndcurves/curve_conversion.h"
#include "ndcurves/exact_cubic.h"
#include "ndcurves/fwd.h"
#include "ndcurves/helpers/effector_spline.h"
#include "ndcurves/helpers/effector_spline_rotation.h"
#include "ndcurves/optimization/definitions.h"
#include "ndcurves/piecewise_curve.h"
#include "ndcurves/polynomial.h"
#include "ndcurves/se3_curve.h"
#include "ndcurves/serialization/curves.hpp"
#include "ndcurves/so3_linear.h"

using namespace std;

namespace ndcurves {
typedef exact_cubic<double, double, true, Eigen::Matrix<double, 1, 1> >
    exact_cubic_one;
typedef exact_cubic_t::spline_constraints spline_constraints_t;

typedef std::pair<double, pointX_t> Waypoint;
typedef std::vector<Waypoint> T_Waypoint;
typedef Eigen::Matrix<double, 1, 1> point_one;
typedef std::pair<double, point_one> WaypointOne;
typedef std::vector<WaypointOne> T_WaypointOne;
typedef std::pair<pointX_t, pointX_t> pair_point_tangent_t;
typedef std::vector<pair_point_tangent_t,
                    Eigen::aligned_allocator<pair_point_tangent_t> >
    t_pair_point_tangent_t;

const double margin = 1e-3;
bool QuasiEqual(const double a, const double b) {
  return std::fabs(a - b) < margin;
}
bool QuasiEqual(const point3_t a, const point3_t b) {
  bool equal = true;
  for (size_t i = 0; i < 3; ++i) {
    equal = equal && QuasiEqual(a[i], b[i]);
  }
  return equal;
}
}  // End namespace ndcurves

using namespace ndcurves;

ostream& operator<<(ostream& os, const point3_t& pt) {
  os << "(" << pt.x() << ", " << pt.y() << ", " << pt.z() << ")";
  return os;
}

void ComparePoints(const transform_t& pt1, const transform_t& pt2,
                   const std::string& errmsg, bool& error,
                   double prec = Eigen::NumTraits<double>::dummy_precision(),
                   bool notequal = false) {
  if (!pt1.isApprox(pt2, prec) && !notequal) {
    error = true;
    std::cout << errmsg << " translation :" << pt1.translation() << " ; "
              << pt2.translation() << std::endl
              << "rotation : " << pt1.rotation() << " ; " << pt2.rotation()
              << std::endl;
  }
}

void ComparePoints(const Eigen::MatrixXd& pt1, const Eigen::MatrixXd& pt2,
                   const std::string& errmsg, bool& error,
                   double prec = Eigen::NumTraits<double>::dummy_precision(),
                   bool notequal = false) {
  if (!pt1.isApprox(pt2, prec) && !(pt1.isZero(prec) && pt2.isZero(prec)) &&
      !notequal) {
    error = true;
    std::cout << errmsg << pt1 << " ; " << pt2 << std::endl;
  }
}

template <typename curve1, typename curve2>
void CompareCurves(const curve1& c1, const curve2& c2,
                   const std::string& errMsg, bool& error,
                   double prec = Eigen::NumTraits<double>::dummy_precision()) {
  double T_min = c1.min();
  double T_max = c1.max();
  if (!QuasiEqual(T_min, c2.min()) || !QuasiEqual(T_max, c2.max())) {
    std::cout
        << errMsg
        << "CompareCurves, ERROR, time min and max of curves do not match ["
        << T_min << "," << T_max << "] "
        << " and [" << c2.min() << "," << c2.max() << "] " << std::endl;
    error = true;
  } else {
    // derivative in T_min and T_max
    ComparePoints(c1.derivate(T_min, 1), c2.derivate(T_min, 1),
                  errMsg + " Derivates at tMin do not match.", error, prec,
                  false);
    ComparePoints(c1.derivate(T_max, 1), c2.derivate(T_max, 1),
                  errMsg + " Derivates at tMax do not match.", error, prec,
                  false);
    // Test values on curves
    for (double i = T_min; i <= T_max; i += 0.01) {
      ComparePoints(c1(i), c2(i),
                    errMsg + " Curves evaluation do not match at t = " +
                        boost::lexical_cast<std::string>(i),
                    error, prec, false);
    }
  }
}

/*Cubic Function tests*/
void PolynomialCubicFunctionTest(bool& error) {
  std::string errMsg("In test CubicFunctionTest ; unexpected result for x ");
  point3_t a(1, 2, 3);
  point3_t b(2, 3, 4);
  point3_t c(3, 4, 5);
  point3_t d(3, 6, 7);
  t_pointX_t vec;
  vec.push_back(a);
  vec.push_back(b);
  vec.push_back(c);
  vec.push_back(d);
  polynomial_t cf(vec.begin(), vec.end(), 0, 1);
  point3_t res1;
  res1 = cf(0);
  point3_t x0(1, 2, 3);
  ComparePoints(x0, res1, errMsg + "(0) ", error);
  point3_t x1(9, 15, 19);
  res1 = cf(1);
  ComparePoints(x1, res1, errMsg + "(1) ", error);
  point3_t x2(3.125, 5.25, 7.125);
  res1 = cf(0.5);
  ComparePoints(x2, res1, errMsg + "(0.5) ", error);
  vec.clear();
  vec.push_back(a);
  vec.push_back(b);
  vec.push_back(c);
  vec.push_back(d);
  polynomial_t cf2(vec, 0.5, 1);
  res1 = cf2(0.5);
  ComparePoints(x0, res1, errMsg + "x3 ", error);
  error = true;
  try {
    cf2(0.4);
  } catch (...) {
    error = false;
  }
  if (error) {
    std::cout << "Evaluation of cubic cf2 error, 0.4 should be an out of range "
                 "value\n";
  }
  error = true;
  try {
    cf2(1.1);
  } catch (...) {
    error = false;
  }
  if (error) {
    std::cout << "Evaluation of cubic cf2 error, 1.1 should be an out of range "
                 "value\n";
  }
  if (!QuasiEqual(cf.max(), 1.0)) {
    error = true;
    std::cout
        << "Evaluation of cubic cf error, MaxBound should be equal to 1\n";
  }
  if (!QuasiEqual(cf.min(), 0.0)) {
    error = true;
    std::cout
        << "Evaluation of cubic cf error, MinBound should be equal to 1\n";
  }
  // Test derivate and compute_derivative
  // Order 1
  curve_abc_t* cf_derivated = cf.compute_derivate_ptr(1);
  ComparePoints(cf.derivate(0, 1), (*cf_derivated)(0),
                errMsg + " - derivate order 1 : ", error);
  ComparePoints(cf.derivate(0.3, 1), (*cf_derivated)(0.3),
                errMsg + " - derivate order 1 : ", error);
  ComparePoints(cf.derivate(0.5, 1), (*cf_derivated)(0.5),
                errMsg + " - derivate order 1 : ", error);
  ComparePoints(cf.derivate(1, 1), (*cf_derivated)(1),
                errMsg + " - derivate order 1 : ", error);
  // Order 2
  polynomial_t cf_derivated_2 = cf.compute_derivate(2);
  ComparePoints(cf.derivate(0, 2), (cf_derivated_2)(0),
                errMsg + " - derivate order 1 : ", error);
  ComparePoints(cf.derivate(0.3, 2), (cf_derivated_2)(0.3),
                errMsg + " - derivate order 1 : ", error);
  ComparePoints(cf.derivate(0.5, 2), (cf_derivated_2)(0.5),
                errMsg + " - derivate order 1 : ", error);
  ComparePoints(cf.derivate(1, 2), (cf_derivated_2)(1),
                errMsg + " - derivate order 1 : ", error);
}

/*bezier_curve Function tests*/
void BezierCurveTest(bool& error) {
  std::string errMsg("In test BezierCurveTest ; unexpected result for x ");
  point3_t a(1, 2, 3);
  point3_t b(2, 3, 4);
  point3_t c(3, 4, 5);
  point3_t d(3, 6, 7);
  std::vector<point3_t> params;
  params.push_back(a);
  // 1d curve in [0,1]
  bezier_t cf1(params.begin(), params.end());
  point3_t res1;
  res1 = cf1(0);
  point3_t x10 = a;
  ComparePoints(x10, res1, errMsg + "1(0) ", error);
  res1 = cf1(1);
  ComparePoints(x10, res1, errMsg + "1(1) ", error);
  // 2d curve in [0,1]
  params.push_back(b);
  bezier_t cf(params.begin(), params.end());
  res1 = cf(0);
  point3_t x20 = a;
  ComparePoints(x20, res1, errMsg + "2(0) ", error);
  point3_t x21 = b;
  res1 = cf(1);
  ComparePoints(x21, res1, errMsg + "2(1) ", error);
  // 3d curve in [0,1]
  params.push_back(c);
  bezier_t cf3(params.begin(), params.end());
  res1 = cf3(0);
  ComparePoints(a, res1, errMsg + "3(0) ", error);
  res1 = cf3(1);
  ComparePoints(c, res1, errMsg + "3(1) ", error);
  // 4d curve in [1,2]
  params.push_back(d);
  bezier_t cf4(params.begin(), params.end(), 1., 2.);
  // testing bernstein polynomials
  bezier_t cf5(params.begin(), params.end(), 1., 2.);
  std::string errMsg2(
      "In test BezierCurveTest ; Bernstein polynomials do not evaluate as "
      "analytical evaluation");
  bezier_t cf5_derivated = cf5.compute_derivate(1);

  for (double d = 1.; d < 2.; d += 0.1) {
    ComparePoints(cf5.evalBernstein(d), cf5(d), errMsg2, error);
    ComparePoints(cf5.evalHorner(d), cf5(d), errMsg2, error);
    ComparePoints(cf5_derivated.evalBernstein(d), cf5_derivated(d), errMsg2,
                  error);
    ComparePoints(cf5_derivated.evalHorner(d), cf5_derivated(d), errMsg2,
                  error);
    ComparePoints(cf5.derivate(d, 1), cf5_derivated(d), errMsg2, error);
  }
  bool error_in(true);
  try {
    cf(-0.4);
  } catch (...) {
    error_in = false;
  }
  if (error_in) {
    std::cout << "Evaluation of bezier cf error, -0.4 should be an out of "
                 "range value\n";
    error = true;
  }
  error_in = true;
  try {
    cf(1.1);
  } catch (...) {
    error_in = false;
  }
  if (error_in) {
    std::cout << "Evaluation of bezier cf error, 1.1 should be an out of range "
                 "value\n";
    error = true;
  }
  if (!QuasiEqual(cf.max(), 1.0)) {
    error = true;
    std::cout
        << "Evaluation of bezier cf error, MaxBound should be equal to 1\n";
  }
  if (!QuasiEqual(cf.min(), 0.0)) {
    error = true;
    std::cout
        << "Evaluation of bezier cf error, MinBound should be equal to 1\n";
  }
}

void BezierCurveTestCompareHornerAndBernstein(bool&)  // error
{
  using namespace std;
  std::vector<double> values;
  for (int i = 0; i < 100000; ++i) {
    values.push_back(rand() / (double)RAND_MAX);
  }
  // first compare regular evaluation (low dim pol)
  point3_t a(1, 2, 3);
  point3_t b(2, 3, 4);
  point3_t c(3, 4, 5);
  point3_t d(3, 6, 7);
  point3_t e(3, 61, 7);
  point3_t f(3, 56, 7);
  point3_t g(3, 36, 7);
  point3_t h(43, 6, 7);
  point3_t i(3, 6, 77);
  std::vector<point3_t> params;
  params.push_back(a);
  params.push_back(b);
  params.push_back(c);
  // 3d curve
  bezier_t cf(params.begin(), params.end());  // defined in [0,1]
  // Check all evaluation of bezier curve
  clock_t s0, e0, s1, e1, s2, e2, s3, e3;
  s0 = clock();
  for (std::vector<double>::const_iterator cit = values.begin();
       cit != values.end(); ++cit) {
    cf(*cit);
  }
  e0 = clock();
  s1 = clock();
  for (std::vector<double>::const_iterator cit = values.begin();
       cit != values.end(); ++cit) {
    cf.evalBernstein(*cit);
  }
  e1 = clock();

  s2 = clock();
  for (std::vector<double>::const_iterator cit = values.begin();
       cit != values.end(); ++cit) {
    cf.evalHorner(*cit);
  }
  e2 = clock();
  s3 = clock();
  for (std::vector<double>::const_iterator cit = values.begin();
       cit != values.end(); ++cit) {
    cf.evalDeCasteljau(*cit);
  }
  e3 = clock();
  std::cout << "time for analytical eval  " << double(e0 - s0) / CLOCKS_PER_SEC
            << std::endl;
  std::cout << "time for bernstein eval   " << double(e1 - s1) / CLOCKS_PER_SEC
            << std::endl;
  std::cout << "time for horner eval      " << double(e2 - s2) / CLOCKS_PER_SEC
            << std::endl;
  std::cout << "time for deCasteljau eval " << double(e3 - s3) / CLOCKS_PER_SEC
            << std::endl;
  std::cout << "now with high order polynomial " << std::endl;
  params.push_back(d);
  params.push_back(e);
  params.push_back(f);
  params.push_back(g);
  params.push_back(h);
  params.push_back(i);
  bezier_t cf2(params.begin(), params.end());
  s1 = clock();
  for (std::vector<double>::const_iterator cit = values.begin();
       cit != values.end(); ++cit) {
    cf2.evalBernstein(*cit);
  }
  e1 = clock();
  s2 = clock();
  for (std::vector<double>::const_iterator cit = values.begin();
       cit != values.end(); ++cit) {
    cf2.evalHorner(*cit);
  }
  e2 = clock();
  s0 = clock();
  for (std::vector<double>::const_iterator cit = values.begin();
       cit != values.end(); ++cit) {
    cf2(*cit);
  }
  e0 = clock();
  s3 = clock();
  for (std::vector<double>::const_iterator cit = values.begin();
       cit != values.end(); ++cit) {
    cf2.evalDeCasteljau(*cit);
  }
  e3 = clock();
  std::cout << "time for analytical eval  " << double(e0 - s0) / CLOCKS_PER_SEC
            << std::endl;
  std::cout << "time for bernstein eval   " << double(e1 - s1) / CLOCKS_PER_SEC
            << std::endl;
  std::cout << "time for horner eval      " << double(e2 - s2) / CLOCKS_PER_SEC
            << std::endl;
  std::cout << "time for deCasteljau eval " << double(e3 - s3) / CLOCKS_PER_SEC
            << std::endl;
}

void BezierDerivativeCurveTest(bool& error) {
  std::string errMsg(
      "In test BezierDerivativeCurveTest ;, Error While checking value of "
      "point on curve : ");
  point3_t a(1, 2, 3);
  point3_t b(2, 3, 4);
  point3_t c(3, 4, 5);
  std::vector<point3_t> params;
  params.push_back(a);
  params.push_back(b);
  params.push_back(c);
  bezier_t cf3(params.begin(), params.end());
  ComparePoints(cf3(0), cf3.derivate(0., 1), errMsg, error, true);
  ComparePoints(point3_t::Zero(), cf3.derivate(0., 100), errMsg, error);
}

void BezierDerivativeCurveTimeReparametrizationTest(bool& error) {
  std::string errMsg(
      "In test BezierDerivativeCurveTimeReparametrizationTest, Error While "
      "checking value of point on curve : ");
  point3_t a(1, 2, 3);
  point3_t b(2, 3, 4);
  point3_t c(3, 4, 5);
  point3_t d(3, 4, 5);
  point3_t e(3, 4, 5);
  point3_t f(3, 4, 5);
  std::vector<point3_t> params;
  params.push_back(a);
  params.push_back(b);
  params.push_back(c);
  params.push_back(d);
  params.push_back(e);
  params.push_back(f);
  double Tmin = 0.;
  double Tmax = 2.;
  double diffT = Tmax - Tmin;
  bezier_t cf(params.begin(), params.end());
  bezier_t cfT(params.begin(), params.end(), Tmin, Tmax);
  ComparePoints(cf(0.5), cfT(1), errMsg, error);
  ComparePoints(cf.derivate(0.5, 1), cfT.derivate(1, 1) * (diffT), errMsg,
                error);
  ComparePoints(cf.derivate(0.5, 2), cfT.derivate(1, 2) * diffT * diffT, errMsg,
                error);
}

void BezierDerivativeCurveConstraintTest(bool& error) {
  std::string errMsg0(
      "In test BezierDerivativeCurveConstraintTest, Error While checking value "
      "of point on curve : ");
  point3_t a(1, 2, 3);
  point3_t b(2, 3, 4);
  point3_t c(3, 4, 5);
  bezier_t::curve_constraints_t constraints(3);
  constraints.init_vel = point3_t(-1, -1, -1);
  constraints.init_acc = point3_t(-2, -2, -2);
  constraints.end_vel = point3_t(-10, -10, -10);
  constraints.end_acc = point3_t(-20, -20, -20);
  std::vector<point3_t> params;
  params.push_back(a);
  params.push_back(b);
  params.push_back(c);
  bezier_t::num_t T_min = 1.0;
  bezier_t::num_t T_max = 3.0;
  bezier_t cf(params.begin(), params.end(), constraints, T_min, T_max);
  ComparePoints(a, cf(T_min), errMsg0, error);
  ComparePoints(c, cf(T_max), errMsg0, error);
  ComparePoints(constraints.init_vel, cf.derivate(T_min, 1), errMsg0, error);
  ComparePoints(constraints.end_vel, cf.derivate(T_max, 1), errMsg0, error);
  ComparePoints(constraints.init_acc, cf.derivate(T_min, 2), errMsg0, error);
  ComparePoints(constraints.end_vel, cf.derivate(T_max, 1), errMsg0, error);
  ComparePoints(constraints.end_acc, cf.derivate(T_max, 2), errMsg0, error);
  std::string errMsg1(
      "In test BezierDerivativeCurveConstraintTest, Error While checking "
      "checking degree of bezier curve :");
  std::string errMsg2(
      "In test BezierDerivativeCurveConstraintTest, Error While checking "
      "checking size of bezier curve :");
  if (cf.degree_ != params.size() + 3) {
    error = true;
    std::cout << errMsg1 << cf.degree_ << " ; " << params.size() + 3
              << std::endl;
  }
  if (cf.size_ != params.size() + 4) {
    error = true;
    std::cout << errMsg2 << cf.size_ << " ; " << params.size() + 4 << std::endl;
  }
}

void toPolynomialConversionTest(bool& error) {
  // bezier to polynomial
  std::string errMsg(
      "In test BezierToPolynomialConversionTest, Error While checking value of "
      "point on curve : ");
  point3_t a(1, 2, 3);
  point3_t b(2, 3, 4);
  point3_t c(3, 4, 5);
  point3_t d(3, 6, 7);
  point3_t e(3, 61, 7);
  point3_t f(3, 56, 7);
  point3_t g(3, 36, 7);
  point3_t h(43, 6, 7);
  point3_t i(3, 6, 77);
  std::vector<point3_t> control_points;
  control_points.push_back(a);
  control_points.push_back(b);
  control_points.push_back(c);
  control_points.push_back(d);
  control_points.push_back(e);
  control_points.push_back(f);
  control_points.push_back(g);
  control_points.push_back(h);
  control_points.push_back(i);
  bezier_t::num_t T_min = 1.0;
  bezier_t::num_t T_max = 3.0;
  bezier_t bc(control_points.begin(), control_points.end(), T_min, T_max);
  polynomial_t pol = polynomial_from_curve<polynomial_t>(bc);
  CompareCurves<polynomial_t, bezier_t>(pol, bc, errMsg, error);
}

void cubicConversionTest(bool& error) {
  std::string errMsg0(
      "In test CubicConversionTest - convert hermite to, Error While checking "
      "value of point on curve : ");
  std::string errMsg1(
      "In test CubicConversionTest - convert bezier to, Error While checking "
      "value of point on curve : ");
  std::string errMsg2(
      "In test CubicConversionTest - convert polynomial to, Error While "
      "checking value of point on curve : ");
  // Create cubic hermite spline : Test hermite to bezier/polynomial
  point3_t p0(1, 2, 3);
  point3_t m0(2, 3, 4);
  point3_t p1(3, 4, 5);
  point3_t m1(3, 6, 7);
  pair_point_tangent_t pair0(p0, m0);
  pair_point_tangent_t pair1(p1, m1);
  t_pair_point_tangent_t control_points;
  control_points.push_back(pair0);
  control_points.push_back(pair1);
  std::vector<double> time_control_points;
  polynomial_t::num_t T_min = 1.0;
  polynomial_t::num_t T_max = 3.0;
  time_control_points.push_back(T_min);
  time_control_points.push_back(T_max);
  cubic_hermite_spline_t chs0(control_points.begin(), control_points.end(),
                              time_control_points);
  // hermite to bezier
  // std::cout<<"======================= \n";
  // std::cout<<"hermite to bezier \n";
  bezier_t bc0 = bezier_from_curve<bezier_t>(chs0);
  CompareCurves<cubic_hermite_spline_t, bezier_t>(chs0, bc0, errMsg0, error);
  // hermite to pol
  // std::cout<<"======================= \n";
  // std::cout<<"hermite to polynomial \n";
  polynomial_t pol0 = polynomial_from_curve<polynomial_t>(chs0);
  CompareCurves<cubic_hermite_spline_t, polynomial_t>(chs0, pol0, errMsg0,
                                                      error);
  // pol to hermite
  // std::cout<<"======================= \n";
  // std::cout<<"polynomial to hermite \n";
  cubic_hermite_spline_t chs1 =
      hermite_from_curve<cubic_hermite_spline_t>(pol0);
  CompareCurves<polynomial_t, cubic_hermite_spline_t>(pol0, chs1, errMsg2,
                                                      error);
  // pol to bezier
  // std::cout<<"======================= \n";
  // std::cout<<"polynomial to bezier \n";
  bezier_t bc1 = bezier_from_curve<bezier_t>(pol0);
  CompareCurves<bezier_t, polynomial_t>(bc1, pol0, errMsg2, error);
  // Bezier to pol
  // std::cout<<"======================= \n";
  // std::cout<<"bezier to polynomial \n";
  polynomial_t pol1 = polynomial_from_curve<polynomial_t>(bc0);
  CompareCurves<bezier_t, polynomial_t>(bc0, pol1, errMsg1, error);
  // bezier => hermite
  // std::cout<<"======================= \n";
  // std::cout<<"bezier to hermite \n";
  cubic_hermite_spline_t chs2 = hermite_from_curve<cubic_hermite_spline_t>(bc0);
  CompareCurves<bezier_t, cubic_hermite_spline_t>(bc0, chs2, errMsg1, error);

  // Test : compute derivative of bezier => Convert it to polynomial
  curve_abc_t* bc_der = bc0.compute_derivate_ptr(1);
  polynomial_t pol_test = polynomial_from_curve<polynomial_t>(*bc_der);
  CompareCurves<curve_abc_t, polynomial_t>(*bc_der, pol_test, errMsg1, error);

  // check that an error is correctly raised when degree > 3:
  point3_t a(1, 2, 3);
  point3_t b(2, 3, 4);
  point3_t c(3, 4, 5);
  point3_t d(3, 6, 7);
  point3_t e(3, 6, 7);
  t_pointX_t vec;
  vec.push_back(a);
  vec.push_back(b);
  vec.push_back(c);
  vec.push_back(d);
  vec.push_back(e);
  polynomial_t pol_4(vec.begin(), vec.end(), 0, 1);
  if (pol_4.degree() != 4) {
    std::cout << "In test CubicConversionTest - Error in the creatin of the "
                 "polynomial"
              << std::endl;
    error = true;
  }
  try {
    cubic_hermite_spline_t chs3 =
        hermite_from_curve<cubic_hermite_spline_t>(pol_4);
    std::cout << "In test CubicConversionTest - Cannot convert to hermite from "
                 "degree > 3, should raise an error"
              << std::endl;
    error = true;
  } catch (std::invalid_argument& /*e*/) {
  }
  try {
    bezier_t b3 = bezier_from_curve<bezier_t>(pol_4);
    std::cout << "In test CubicConversionTest - Cannot convert to bezier from "
                 "degree > 3, should raise an error"
              << std::endl;
    error = true;
  } catch (std::invalid_argument& /*e*/) {
  }
}

/*Exact Cubic Function tests*/
void ExactCubicNoErrorTest(bool& error) {
  // Create an exact cubic spline with 7 waypoints => 6 polynomials defined in
  // [0.0,3.0]
  ndcurves::T_Waypoint waypoints;
  for (double i = 0.0; i <= 3.0; i = i + 0.5) {
    waypoints.push_back(std::make_pair(i, point3_t(i, i, i)));
  }
  exact_cubic_t exactCubic(waypoints.begin(), waypoints.end());
  // Test number of polynomials in exact cubic
  std::size_t numberSegments = exactCubic.getNumberSplines();
  if (numberSegments != 6) {
    error = true;
    std::cout
        << "In ExactCubicNoErrorTest, Error While checking number of splines"
        << numberSegments << " ; " << 6 << std::endl;
  }
  // Test getSplineAt function
  for (std::size_t i = 0; i < numberSegments; i++) {
    exactCubic.getSplineAt(i);
  }
  // Other tests
  try {
    exactCubic(0.0);
    exactCubic(3.0);
  } catch (...) {
    error = true;
    std::cout << "Evaluation of ExactCubicNoErrorTest error when testing value "
                 "on bounds\n";
  }
  error = true;
  try {
    exactCubic(3.2);
  } catch (...) {
    error = false;
  }
  if (error) {
    std::cout << "Evaluation of exactCubic cf error, 3.2 should be an out of "
                 "range value\n";
  }
  if (!QuasiEqual(exactCubic.max(), 3.0)) {
    error = true;
    std::cout << "Evaluation of exactCubic error, MaxBound should be equal to "
                 "3 but is : "
              << exactCubic.max() << "\n";
  }
  if (!QuasiEqual(exactCubic.min(), 0.0)) {
    error = true;
    std::cout << "Evaluation of exactCubic error, MinBound should be equal to "
                 "0 but is : "
              << exactCubic.min() << "\n";
  }
}

/*Exact Cubic Function tests*/
void ExactCubicTwoPointsTest(bool& error) {
  // Create an exact cubic spline with 2 waypoints => 1 polynomial defined in
  // [0.0,1.0]
  ndcurves::T_Waypoint waypoints;
  for (double i = 0.0; i < 2.0; ++i) {
    waypoints.push_back(std::make_pair(i, point3_t(i, i, i)));
  }
  exact_cubic_t exactCubic(waypoints.begin(), waypoints.end());
  point3_t res1 = exactCubic(0);
  std::string errmsg0(
      "in ExactCubicTwoPointsTest, Error While checking that given wayPoints  "
      "are crossed (expected / obtained)");
  ComparePoints(point3_t(0, 0, 0), res1, errmsg0, error);
  res1 = exactCubic(1);
  ComparePoints(point3_t(1, 1, 1), res1, errmsg0, error);
  // Test number of polynomials in exact cubic
  std::size_t numberSegments = exactCubic.getNumberSplines();
  if (numberSegments != 1) {
    error = true;
    std::cout
        << "In ExactCubicTwoPointsTest, Error While checking number of splines"
        << numberSegments << " ; " << 1 << std::endl;
  }
  // Test getSplineAt
  std::string errmsg1(
      "in ExactCubicTwoPointsTest, Error While checking value on curve");
  ComparePoints(exactCubic(0.5), (exactCubic.getSplineAt(0))(0.5), errmsg1,
                error);
}

void ExactCubicOneDimTest(bool& error) {
  ndcurves::T_WaypointOne waypoints;
  point_one zero;
  zero(0, 0) = 9;
  point_one one;
  one(0, 0) = 14;
  point_one two;
  two(0, 0) = 25;
  waypoints.push_back(std::make_pair(0., zero));
  waypoints.push_back(std::make_pair(1., one));
  waypoints.push_back(std::make_pair(2., two));
  exact_cubic_one exactCubic(waypoints.begin(), waypoints.end());
  point_one res1 = exactCubic(0);
  std::string errmsg(
      "in ExactCubicOneDim Error While checking that given wayPoints  are "
      "crossed (expected / obtained)");
  ComparePoints(zero, res1, errmsg, error);
  res1 = exactCubic(1);
  ComparePoints(one, res1, errmsg, error);
}

void CheckWayPointConstraint(
    const std::string& errmsg, const double step, const ndcurves::T_Waypoint&,
    const exact_cubic_t* curve, bool& error,
    double prec = Eigen::NumTraits<double>::dummy_precision()) {
  point3_t res1;
  for (double i = 0; i <= 1; i = i + step) {
    res1 = (*curve)(i);
    ComparePoints(point3_t(i, i, i), res1, errmsg, error, prec);
  }
}

void ExactCubicPointsCrossedTest(bool& error) {
  ndcurves::T_Waypoint waypoints;
  for (double i = 0; i <= 1; i = i + 0.2) {
    waypoints.push_back(std::make_pair(i, point3_t(i, i, i)));
  }
  exact_cubic_t exactCubic(waypoints.begin(), waypoints.end());
  std::string errmsg(
      "Error While checking that given wayPoints are crossed (expected / "
      "obtained)");
  CheckWayPointConstraint(errmsg, 0.2, waypoints, &exactCubic, error);
}

void ExactCubicVelocityConstraintsTest(bool& error) {
  ndcurves::T_Waypoint waypoints;
  for (double i = 0; i <= 1; i = i + 0.2) {
    waypoints.push_back(std::make_pair(i, point3_t(i, i, i)));
  }
  std::string errmsg(
      "Error in ExactCubicVelocityConstraintsTest (1); while checking that "
      "given wayPoints are crossed (expected / "
      "obtained)");
  spline_constraints_t constraints(3);
  constraints.end_vel = point3_t(0, 0, 0);
  constraints.init_vel = point3_t(0, 0, 0);
  constraints.end_acc = point3_t(0, 0, 0);
  constraints.init_acc = point3_t(0, 0, 0);
  exact_cubic_t exactCubic(waypoints.begin(), waypoints.end(), constraints);
  // now check that init and end velocity are 0
  CheckWayPointConstraint(errmsg, 0.2, waypoints, &exactCubic, error);
  std::string errmsg3(
      "Error in ExactCubicVelocityConstraintsTest (2); while checking "
      "derivative (expected / obtained)");
  // now check derivatives
  ComparePoints(constraints.init_vel, exactCubic.derivate(0, 1), errmsg3, error,
                1e-10);
  ComparePoints(constraints.end_vel, exactCubic.derivate(1, 1), errmsg3, error,
                1e-10);
  ComparePoints(constraints.init_acc, exactCubic.derivate(0, 2), errmsg3, error,
                1e-10);
  ComparePoints(constraints.end_acc, exactCubic.derivate(1, 2), errmsg3, error,
                1e-10);
  constraints.end_vel = point3_t(1, 2, 3);
  constraints.init_vel = point3_t(-1, -2, -3);
  constraints.end_acc = point3_t(4, 5, 6);
  constraints.init_acc = point3_t(-4, -4, -6);
  std::string errmsg2(
      "Error in ExactCubicVelocityConstraintsTest (3); while checking that "
      "given wayPoints are crossed (expected / "
      "obtained)");
  exact_cubic_t exactCubic2(waypoints.begin(), waypoints.end(), constraints);
  CheckWayPointConstraint(errmsg2, 0.2, waypoints, &exactCubic2, error, 1e-10);
  std::string errmsg4(
      "Error in ExactCubicVelocityConstraintsTest (4); while checking "
      "derivative (expected / obtained)");
  // now check derivatives
  ComparePoints(constraints.init_vel, exactCubic2.derivate(0, 1), errmsg4,
                error, 1e-10);
  ComparePoints(constraints.end_vel, exactCubic2.derivate(1, 1), errmsg4, error,
                1e-10);
  ComparePoints(constraints.init_acc, exactCubic2.derivate(0, 2), errmsg4,
                error, 1e-10);
  ComparePoints(constraints.end_acc, exactCubic2.derivate(1, 2), errmsg4, error,
                1e-10);
}

template <typename CurveType>
void CheckPointOnline(const std::string& errmsg, const point3_t& A,
                      const point3_t& B, const double target,
                      const CurveType* curve, bool& error) {
  point3_t res1 = curve->operator()(target);
  point3_t ar = (res1 - A);
  ar.normalize();
  point3_t rb = (B - res1);
  rb.normalize();
  if (ar.dot(rb) < 0.99999) {
    error = true;
    std::cout << errmsg << " ; " << A.transpose() << "\n ; " << B.transpose()
              << "\n ; " << target << " ; " << res1.transpose() << std::endl;
  }
}

void EffectorTrajectoryTest(bool& error) {
  // create arbitrary trajectory
  ndcurves::T_Waypoint waypoints;
  for (double i = 0; i <= 10; i = i + 2) {
    waypoints.push_back(std::make_pair(i, point3_t(i, i, i)));
  }
  helpers::exact_cubic_t* eff_traj = helpers::effector_spline(
      waypoints.begin(), waypoints.end(), Eigen::Vector3d::UnitZ(),
      Eigen::Vector3d(0, 0, 2), 1, 0.02, 1, 0.5);
  point3_t zero(0, 0, 0);
  point3_t off1(0, 0, 1);
  point3_t off2(10, 10, 10.02);
  point3_t end(10, 10, 10);
  std::string errmsg(
      "Error in EffectorTrajectoryTest; while checking waypoints (expected / "
      "obtained)");
  std::string errmsg2(
      "Error in EffectorTrajectoryTest; while checking derivative (expected / "
      "obtained)");
  // first check start / goal positions
  ComparePoints(zero, (*eff_traj)(0), errmsg, error);
  ComparePoints(off1, (*eff_traj)(1), errmsg, error);
  ComparePoints(off2, (*eff_traj)(9.5), errmsg, error);
  ComparePoints(end, (*eff_traj)(10), errmsg, error);
  // now check derivatives
  ComparePoints(zero, (*eff_traj).derivate(0, 1), errmsg2, error);
  ComparePoints(zero, (*eff_traj).derivate(10, 1), errmsg2, error);
  ComparePoints(zero, (*eff_traj).derivate(0, 2), errmsg2, error);
  ComparePoints(zero, (*eff_traj).derivate(10, 2), errmsg2, error);
  // check that end and init splines are line
  std::string errmsg3(
      "Error in EffectorTrajectoryTest; while checking that init/end splines "
      "are line (point A/ point B, time value / "
      "point obtained) \n");
  for (double i = 0.1; i < 1; i += 0.1) {
    CheckPointOnline<helpers::exact_cubic_t>(
        errmsg3, (*eff_traj)(0), (*eff_traj)(1), i, eff_traj, error);
  }
  for (double i = 9.981; i < 10; i += 0.002) {
    CheckPointOnline<helpers::exact_cubic_t>(
        errmsg3, (*eff_traj)(9.5), (*eff_traj)(10), i, eff_traj, error);
  }
  delete eff_traj;
}

helpers::quat_t GetXRotQuat(const double theta) {
  Eigen::AngleAxisd m(theta, Eigen::Vector3d::UnitX());
  return helpers::quat_t(Eigen::Quaterniond(m).coeffs().data());
}

double GetXRotFromQuat(helpers::quat_ref_const_t q) {
  Eigen::Quaterniond quat(q.data());
  Eigen::AngleAxisd m(quat);
  return m.angle() / M_PI * 180.;
}

void EffectorSplineRotationNoRotationTest(bool& error) {
  // create arbitrary trajectory
  ndcurves::T_Waypoint waypoints;
  for (double i = 0; i <= 10; i = i + 2) {
    waypoints.push_back(std::make_pair(i, point3_t(i, i, i)));
  }
  helpers::effector_spline_rotation eff_traj(waypoints.begin(),
                                             waypoints.end());
  helpers::config_t q_init;
  q_init << 0., 0., 0., 0., 0., 0., 1.;
  helpers::config_t q_end;
  q_end << 10., 10., 10., 0., 0., 0., 1.;
  helpers::config_t q_to;
  q_to << 0., 0, 0.02, 0., 0., 0., 1.;
  helpers::config_t q_land;
  q_land << 10, 10, 10.02, 0, 0., 0., 1.;
  helpers::config_t q_mod;
  q_mod << 6., 6., 6., 0., 0., 0., 1.;
  std::string errmsg(
      "Error in EffectorSplineRotationNoRotationTest; while checking waypoints "
      "(expected / obtained)");
  ComparePoints(q_init, eff_traj(0), errmsg, error);
  ComparePoints(q_to, eff_traj(0.02), errmsg, error);
  ComparePoints(q_land, eff_traj(9.98), errmsg, error);
  ComparePoints(q_mod, eff_traj(6), errmsg, error);
  ComparePoints(q_end, eff_traj(10), errmsg, error);
}

void EffectorSplineRotationRotationTest(bool& error) {
  // create arbitrary trajectory
  ndcurves::T_Waypoint waypoints;
  for (double i = 0; i <= 10; i = i + 2) {
    waypoints.push_back(std::make_pair(i, point3_t(i, i, i)));
  }
  helpers::quat_t init_quat = GetXRotQuat(M_PI);
  helpers::effector_spline_rotation eff_traj(waypoints.begin(), waypoints.end(),
                                             init_quat);
  helpers::config_t q_init = helpers::config_t::Zero();
  q_init.tail<4>() = init_quat;
  helpers::config_t q_end;
  q_end << 10., 10., 10., 0., 0., 0., 1.;
  helpers::config_t q_to = q_init;
  q_to(2) += 0.02;
  helpers::config_t q_land = q_end;
  q_land(2) += 0.02;
  helpers::quat_t q_mod = GetXRotQuat(M_PI_2);
  ;
  std::string errmsg(
      "Error in EffectorSplineRotationRotationTest; while checking waypoints "
      "(expected / obtained)");
  ComparePoints(q_init, eff_traj(0), errmsg, error);
  ComparePoints(q_to, eff_traj(0.02), errmsg, error);
  ComparePoints(q_land, eff_traj(9.98), errmsg, error);
  ComparePoints(q_mod, eff_traj(5).tail<4>(), errmsg, error);
  ComparePoints(q_end, eff_traj(10), errmsg, error);
}

void EffectorSplineRotationWayPointRotationTest(bool& error) {
  // create arbitrary trajectory
  ndcurves::T_Waypoint waypoints;
  for (double i = 0; i <= 10; i = i + 2) {
    waypoints.push_back(std::make_pair(i, point3_t(i, i, i)));
  }
  helpers::quat_t init_quat = GetXRotQuat(0);
  helpers::t_waypoint_quat_t quat_waypoints_;
  helpers::quat_t q_pi_0 = GetXRotQuat(0);
  helpers::quat_t q_pi_2 = GetXRotQuat(M_PI_2);
  helpers::quat_t q_pi = GetXRotQuat(M_PI);
  quat_waypoints_.push_back(std::make_pair(0.4, q_pi_0));
  quat_waypoints_.push_back(std::make_pair(6, q_pi_2));
  quat_waypoints_.push_back(std::make_pair(8, q_pi));
  helpers::effector_spline_rotation eff_traj(waypoints.begin(), waypoints.end(),
                                             quat_waypoints_.begin(),
                                             quat_waypoints_.end());
  helpers::config_t q_init = helpers::config_t::Zero();
  q_init.tail<4>() = init_quat;
  helpers::config_t q_end;
  q_end << 10., 10., 10., 0., 0., 0., 1.;
  q_end.tail<4>() = q_pi;
  helpers::config_t q_mod;
  q_mod.head<3>() = point3_t(6, 6, 6);
  q_mod.tail<4>() = q_pi_2;
  helpers::config_t q_to = q_init;
  q_to(2) += 0.02;
  helpers::config_t q_land = q_end;
  q_land(2) += 0.02;
  std::string errmsg(
      "Error in EffectorSplineRotationWayPointRotationTest; while checking "
      "waypoints (expected / obtained)");
  ComparePoints(q_init, eff_traj(0), errmsg, error);
  ComparePoints(q_to, eff_traj(0.02), errmsg, error);
  ComparePoints(q_land, eff_traj(9.98), errmsg, error);
  ComparePoints(q_mod, eff_traj(6), errmsg, error);
  ComparePoints(q_end, eff_traj(10), errmsg, error);
}

void TestReparametrization(bool& error) {
  helpers::rotation_spline s;
  const helpers::exact_cubic_constraint_one_dim& sp = s.time_reparam_;
  if (!QuasiEqual(sp.min(), 0.0)) {
    std::cout << "in TestReparametrization; min value is not 0, got "
              << sp.min() << std::endl;
    error = true;
  }
  if (!QuasiEqual(sp.max(), 1.0)) {
    std::cout << "in TestReparametrization; max value is not 1, got "
              << sp.max() << std::endl;
    error = true;
  }
  if (!QuasiEqual(sp(1)[0], 1.0)) {
    std::cout << "in TestReparametrization; end value is not 1, got "
              << sp(1)[0] << std::endl;
    error = true;
  }
  if (!QuasiEqual(sp(0)[0], 0.0)) {
    std::cout << "in TestReparametrization; init value is not 0, got "
              << sp(0)[0] << std::endl;
    error = true;
  }
  for (double i = 0; i < 1; i += 0.002) {
    if (sp(i)[0] > sp(i + 0.002)[0]) {
      std::cout << "in TestReparametrization; reparametrization not monotonous "
                << sp.max() << std::endl;
      error = true;
    }
  }
}

point3_t randomPoint(const double min, const double max) {
  point3_t p;
  for (size_t i = 0; i < 3; ++i) {
    p[i] = (rand() / (double)RAND_MAX) * (max - min) + min;
  }
  return p;
}

void BezierEvalDeCasteljau(bool& error) {
  using namespace std;
  std::vector<double> values;
  for (int i = 0; i < 100000; ++i) {
    values.push_back(rand() / RAND_MAX);
  }
  // first compare regular evaluation (low dim pol)
  point3_t a(1, 2, 3);
  point3_t b(2, 3, 4);
  point3_t c(3, 4, 5);
  point3_t d(3, 6, 7);
  point3_t e(3, 61, 7);
  point3_t f(3, 56, 7);
  point3_t g(3, 36, 7);
  point3_t h(43, 6, 7);
  point3_t i(3, 6, 77);
  std::vector<point3_t> params;
  params.push_back(a);
  params.push_back(b);
  params.push_back(c);
  // 3d curve
  bezier_t cf(params.begin(), params.end());
  std::string errmsg(
      "Error in BezierEvalDeCasteljau; while comparing actual bezier "
      "evaluation and de Casteljau : ");
  for (std::vector<double>::const_iterator cit = values.begin();
       cit != values.end(); ++cit) {
    ComparePoints(cf.evalDeCasteljau(*cit), cf(*cit), errmsg, error);
  }
  params.push_back(d);
  params.push_back(e);
  params.push_back(f);
  params.push_back(g);
  params.push_back(h);
  params.push_back(i);
  bezier_t cf2(params.begin(), params.end());
  for (std::vector<double>::const_iterator cit = values.begin();
       cit != values.end(); ++cit) {
    ComparePoints(cf2.evalDeCasteljau(*cit), cf2(*cit), errmsg, error);
  }
}

void BezierElevate(bool& error) {
  using namespace std;
  std::vector<double> values;
  for (int i = 0; i < 10; ++i) {
    values.push_back(double(rand()) / double(RAND_MAX));
  }
  // first compare regular evaluation (low dim pol)
  point3_t a(1, 2, 3);
  point3_t b(2, 3, 4);
  point3_t c(3, 4, 5);
  std::vector<point3_t> params;
  params.push_back(a);
  params.push_back(b);
  params.push_back(c);
  // 3d curve
  bezier_t cf(params.begin(), params.end());
  bezier_t cf2 = cf.elevate(1);
  bezier_t cf3 = cf2.elevate(1);
  if (cf2.degree() - cf.degree() != 1 && cf3.degree() - cf.degree() != 2) {
    error = true;
    std::string errmsg(
        "Error in BezierElevate; Degree mismatched for elevated curves. "
        "Expected 1 / 2, got:  ");
    std::cout << errmsg << cf2.degree() - cf.degree() << " ; "
              << cf3.degree() - cf.degree() << std::endl;
  }
  std::string errmsg(
      "Error in BezierElevate; Elevated curves do not have the same evaluation "
      ": ");
  for (std::vector<double>::const_iterator cit = values.begin();
       cit != values.end() - 1; ++cit) {
    ComparePoints(cf2(*cit), cf(*(cit)), errmsg, error);
    ComparePoints(cf3(*cit), cf(*cit), errmsg, error);
  }
}

/**
 * @brief BezierSplitCurve test the 'split' method of bezier curve
 * @param error
 */
void BezierSplitCurve(bool& error) {
  // test for degree 5
  size_t n = 5;
  double t_min = 0.2;
  double t_max = 10;
  double aux0, aux1;
  std::string errMsg0(
      "BezierSplitCurve, ERROR initial point of the splitted curve doesn't "
      "correspond to the original");
  std::string errMsg1(
      "BezierSplitCurve, ERROR splitting point of the splitted curve doesn't "
      "correspond to the original");
  std::string errMsg2(
      "BezierSplitCurve, ERROR final point of the splitted curve doesn't "
      "correspond to the original");
  std::string errMsg3(
      "BezierSplitCurve, ERROR while checking value on curve and curves "
      "splitted");
  std::string errMsg4(
      "BezierSplitCurve, ERROR Degree of the splitted curve are not the same "
      "as the original curve");
  std::string errMsg5(
      "BezierSplitCurve, ERROR duration of the splitted curve doesn't "
      "correspond to the original");
  std::string errMsg6(
      "BezierSplitCurve, ERROR while checking value on curve extracted");
  for (size_t i = 0; i < 1; ++i) {
    // build a random curve and split it at random time :
    // std::cout<<"build a random curve"<<std::endl;
    point3_t a;
    std::vector<point3_t> wps;
    for (size_t j = 0; j <= n; ++j) {
      wps.push_back(randomPoint(-10., 10.));
    }
    double t0 = (rand() / (double)RAND_MAX) * (t_max - t_min) + t_min;
    double t1 = (rand() / (double)RAND_MAX) * (t_max - t0) + t0;
    double ts = (rand() / (double)RAND_MAX) * (t1 - t0) + t0;
    bezier_t c(wps.begin(), wps.end(), t0, t1);
    std::pair<bezier_t, bezier_t> cs = c.split(ts);
    // test on splitted curves :
    if (!((c.degree_ == cs.first.degree_) &&
          (c.degree_ == cs.second.degree_))) {
      error = true;
      std::cout << errMsg4 << std::endl;
    }
    aux0 = c.max() - c.min();
    aux1 =
        (cs.first.max() - cs.first.min() + cs.second.max() - cs.second.min());
    if (!QuasiEqual(aux0, aux1)) {
      error = true;
      std::cout << errMsg5 << std::endl;
    }
    if (!QuasiEqual(cs.first.max(), ts)) {
      error = true;
      std::cout << errMsg0 << std::endl;
    }
    ComparePoints(c(t0), cs.first(t0), errMsg0, error);
    ComparePoints(cs.first(ts), cs.second(ts), errMsg1, error);
    ComparePoints(c(t1), cs.second(cs.second.max()), errMsg2, error);
    // check along curve :
    double ti = t0;
    while (ti <= ts) {
      ComparePoints(cs.first(ti), c(ti), errMsg3, error);
      ti += 0.01;
    }
    while (ti <= t1) {
      ComparePoints(cs.second(ti), c(ti), errMsg3, error);
      ti += 0.01;
    }
    // Test extract function
    bezier_t bezier_extracted0 =
        c.extract(t0 + 0.01, t1 - 0.01);  // T_min < t0 < t1 < T_max
    for (double t = bezier_extracted0.min(); t < bezier_extracted0.max();
         t += 0.01) {
      ComparePoints(bezier_extracted0(t), c(t), errMsg6, error);
    }
    bezier_t bezier_extracted1 =
        c.extract(t0, t1 - 0.01);  // T_min = t0 < t1 < T_max
    for (double t = bezier_extracted1.min(); t < bezier_extracted1.max();
         t += 0.01) {
      ComparePoints(bezier_extracted1(t), c(t), errMsg6, error);
    }
    bezier_t bezier_extracted2 =
        c.extract(t0 + 0.01, t1);  // T_min < t0 < t1 = T_max
    for (double t = bezier_extracted2.min(); t < bezier_extracted2.max();
         t += 0.01) {
      ComparePoints(bezier_extracted2(t), c(t), errMsg6, error);
    }
    bezier_t bezier_extracted3 = c.extract(t0, t1);  // T_min = t0 < t1 = T_max
    for (double t = bezier_extracted3.min(); t < bezier_extracted3.max();
         t += 0.01) {
      ComparePoints(bezier_extracted3(t), c(t), errMsg6, error);
    }
  }
}

/* cubic hermite spline function test */
void CubicHermitePairsPositionDerivativeTest(bool& error) {
  try {
    std::string errmsg1(
        "in Cubic Hermite 2 pairs (pos,vel), Error While checking that given "
        "wayPoints are crossed (expected / "
        "obtained) : ");
    std::string errmsg2(
        "in Cubic Hermite 2 points, Error While checking value of point on "
        "curve : ");
    std::string errmsg3(
        "in Cubic Hermite 2 points, Error While checking value of tangent on "
        "curve : ");
    std::vector<pair_point_tangent_t> control_points;
    point3_t res1;
    point3_t p0(0., 0., 0.);
    point3_t p1(1., 2., 3.);
    point3_t p2(4., 4., 4.);
    point3_t t0(0.5, 0.5, 0.5);
    point3_t t1(0.1, 0.2, -0.5);
    point3_t t2(0.1, 0.2, 0.3);
    std::vector<double> time_control_points, time_control_points_test;
    // Two pairs
    control_points.clear();
    control_points.push_back(pair_point_tangent_t(p0, t0));
    control_points.push_back(pair_point_tangent_t(p1, t1));
    time_control_points.push_back(0.);  // Time at P0
    time_control_points.push_back(1.);  // Time at P1
    // Create cubic hermite spline
    cubic_hermite_spline_t cubic_hermite_spline_1Pair(
        control_points.begin(), control_points.end(), time_control_points);
    // Dimension
    if (cubic_hermite_spline_1Pair.dim() != 3) {
      error = true;
      std::cout
          << "Cubic hermite spline test, Error : Dimension of curve is wrong\n";
    }
    // Check
    res1 = cubic_hermite_spline_1Pair(0.);  // t=0
    ComparePoints(p0, res1, errmsg1, error);
    res1 = cubic_hermite_spline_1Pair(1.);  // t=1
    ComparePoints(p1, res1, errmsg1, error);
    // Test derivative : two pairs
    res1 = cubic_hermite_spline_1Pair.derivate(0., 1);
    ComparePoints(t0, res1, errmsg3, error);
    res1 = cubic_hermite_spline_1Pair.derivate(1., 1);
    ComparePoints(t1, res1, errmsg3, error);
    // Three pairs
    control_points.push_back(pair_point_tangent_t(p2, t2));
    time_control_points.clear();
    time_control_points.push_back(0.);  // Time at P0
    time_control_points.push_back(2.);  // Time at P1
    time_control_points.push_back(5.);  // Time at P2
    cubic_hermite_spline_t cubic_hermite_spline_2Pairs(
        control_points.begin(), control_points.end(), time_control_points);
    // Check
    res1 = cubic_hermite_spline_2Pairs(0.);  // t=0
    ComparePoints(p0, res1, errmsg1, error);
    res1 = cubic_hermite_spline_2Pairs(2.);  // t=2
    ComparePoints(p1, res1, errmsg2, error);
    res1 = cubic_hermite_spline_2Pairs(5.);  // t=5
    ComparePoints(p2, res1, errmsg1, error);
    // Test derivative : three pairs
    res1 = cubic_hermite_spline_2Pairs.derivate(0., 1);
    ComparePoints(t0, res1, errmsg3, error);
    res1 = cubic_hermite_spline_2Pairs.derivate(2., 1);
    ComparePoints(t1, res1, errmsg3, error);
    res1 = cubic_hermite_spline_2Pairs.derivate(5., 1);
    ComparePoints(t2, res1, errmsg3, error);
    // Test time control points by default [0,1] => with N control points :
    // Time at P0= 0. | Time at P1= 1.0/(N-1) | Time at P2= 2.0/(N-1) | ... |
    // Time at P_(N-1)= (N-1)/(N-1)= 1.0
    time_control_points_test.clear();
    time_control_points_test.push_back(0.);   // Time at P0
    time_control_points_test.push_back(0.5);  // Time at P1
    time_control_points_test.push_back(1.0);  // Time at P2
    cubic_hermite_spline_2Pairs.setTime(time_control_points_test);
    res1 = cubic_hermite_spline_2Pairs(0.);  // t=0
    ComparePoints(p0, res1, errmsg1, error);
    res1 = cubic_hermite_spline_2Pairs(0.5);  // t=0.5
    ComparePoints(p1, res1, errmsg2, error);
    res1 = cubic_hermite_spline_2Pairs(1.);  // t=1
    ComparePoints(p2, res1, errmsg1, error);
    // Test getTime
    try {
      cubic_hermite_spline_2Pairs.getTime();
    } catch (...) {
      error = false;
    }
    if (error) {
      std::cout << "Cubic hermite spline test, Error when calling getTime\n";
    }
    // Test derivative : three pairs, time default
    res1 = cubic_hermite_spline_2Pairs.derivate(0., 1);
    ComparePoints(t0, res1, errmsg3, error);
    res1 = cubic_hermite_spline_2Pairs.derivate(0.5, 1);
    ComparePoints(t1, res1, errmsg3, error);
    res1 = cubic_hermite_spline_2Pairs.derivate(1., 1);
    ComparePoints(t2, res1, errmsg3, error);

    pointX_t p_derivate, p_compute_derivate;
    for (size_t order = 1; order < 5; ++order) {
      std::stringstream ss;
      ss << "in Cubic Hermite 2 points, "
            "compute_derivate do not lead to the same results as derivate for "
            "order = ";
      ss << order << std::endl;
      curve_ptr_t derivate_ptr(
          cubic_hermite_spline_2Pairs.compute_derivate_ptr(order));
      double t = 0.;
      while (t <= 1.) {
        p_derivate = cubic_hermite_spline_2Pairs.derivate(t, order);
        p_compute_derivate = derivate_ptr->operator()(t);
        ComparePoints(p_derivate, p_compute_derivate, ss.str(), error);
        t += 0.1;
      }
    }
  } catch (...) {
    error = true;
    std::cout << "Error in CubicHermitePairsPositionDerivativeTest"
              << std::endl;
  }
}

void piecewiseCurveTest(bool& error) {
  try {
    // TEST WITH POLYNOMIALS
    std::string errmsg1(
        "in piecewise polynomial curve test, Error While checking value of "
        "point on curve : ");
    point3_t a(1, 1, 1);  // in [0,1[
    point3_t b(2, 1, 1);  // in [1,2[
    point3_t c(3, 1, 1);  // in [2,3]
    point3_t res;
    t_pointX_t vec1, vec2, vec3;
    vec1.push_back(a);  // x=1, y=1, z=1
    vec2.push_back(b);  // x=2, y=1, z=1
    vec3.push_back(c);  // x=3, y=1, z=1
    // Create three polynomials of constant value in the interval of definition
    std::shared_ptr<polynomial_t> pol1_ptr =
        std::make_shared<polynomial_t>(vec1.begin(), vec1.end(), 0, 1);
    std::shared_ptr<polynomial_t> pol2_ptr =
        std::make_shared<polynomial_t>(vec2.begin(), vec2.end(), 1, 2);
    std::shared_ptr<polynomial_t> pol3_ptr =
        std::make_shared<polynomial_t>(vec3.begin(), vec3.end(), 2, 3);
    // 1 polynomial in curve
    piecewise_t pc(pol1_ptr);
    res = pc(0.5);
    ComparePoints(a, res, errmsg1, error);
    // 3 polynomials in curve
    pc.add_curve_ptr(pol2_ptr);
    pc.add_curve_ptr(pol3_ptr);
    // Check values on piecewise curve
    // t in [0,1[ -> res=a
    res = pc(0.);
    ComparePoints(a, res, errmsg1, error);
    res = pc(0.5);
    ComparePoints(a, res, errmsg1, error);
    // t in [1,2[ -> res=b
    res = pc(1.0);
    ComparePoints(b, res, errmsg1, error);
    res = pc(1.5);
    ComparePoints(b, res, errmsg1, error);
    // t in [2,3] -> res=c
    res = pc(2.0);
    ComparePoints(c, res, errmsg1, error);
    res = pc(3.0);
    ComparePoints(c, res, errmsg1, error);
    // Create piecewise curve C0 from bezier
    point3_t a0(1, 2, 3);
    point3_t b0(2, 3, 4);
    point3_t c0(3, 4, 5);
    point3_t d0(4, 5, 6);
    std::vector<point3_t> params0;
    std::vector<point3_t> params1;
    params0.push_back(a0);  // bezier between [0,1]
    params0.push_back(b0);
    params0.push_back(c0);
    params0.push_back(d0);
    params1.push_back(d0);  // bezier between [1,2]
    params1.push_back(c0);
    params1.push_back(b0);
    params1.push_back(a0);
    std::shared_ptr<bezier_t> bc0_ptr =
        std::make_shared<bezier_t>(params0.begin(), params0.end(), 0., 1.);
    std::shared_ptr<bezier_t> bc1_ptr =
        std::make_shared<bezier_t>(params1.begin(), params1.end(), 1., 2.);
    piecewise_t pc_C0(bc0_ptr);
    pc_C0.add_curve_ptr(bc1_ptr);
    // Check value in t=0.5 and t=1.5
    res = pc_C0(0.0);
    ComparePoints(a0, res, errmsg1, error);
    res = pc_C0(1.0);
    ComparePoints(d0, res, errmsg1, error);
    res = pc_C0(2.0);
    ComparePoints(a0, res, errmsg1, error);
    // Create piecewise curve C1 from Hermite
    point3_t p0(0., 0., 0.);
    point3_t p1(1., 2., 3.);
    point3_t p2(4., 4., 4.);
    point3_t t0(0.5, 0.5, 0.5);
    point3_t t1(0.1, 0.2, -0.5);
    point3_t t2(0.1, 0.2, 0.3);
    std::vector<pair_point_tangent_t> control_points_0;
    control_points_0.push_back(pair_point_tangent_t(p0, t0));
    control_points_0.push_back(
        pair_point_tangent_t(p1, t1));  // control_points_0 = 1st piece of curve
    std::vector<pair_point_tangent_t> control_points_1;
    control_points_1.push_back(pair_point_tangent_t(p1, t1));
    control_points_1.push_back(
        pair_point_tangent_t(p2, t2));  // control_points_1 = 2nd piece of curve
    std::vector<double> time_control_points0, time_control_points1;
    time_control_points0.push_back(0.);
    time_control_points0.push_back(1.);  // hermite 0 between [0,1]
    time_control_points1.push_back(1.);
    time_control_points1.push_back(3.);  // hermite 1 between [1,3]
    std::shared_ptr<cubic_hermite_spline_t> chs0_ptr =
        std::make_shared<cubic_hermite_spline_t>(control_points_0.begin(),
                                                 control_points_0.end(),
                                                 time_control_points0);
    std::shared_ptr<cubic_hermite_spline_t> chs1_ptr =
        std::make_shared<cubic_hermite_spline_t>(control_points_1.begin(),
                                                 control_points_1.end(),
                                                 time_control_points1);
    piecewise_t pc_C1(chs0_ptr);
    pc_C1.add_curve_ptr(chs1_ptr);
    // Create piecewise curve C2
    point3_t a1(0, 0, 0);
    point3_t b1(1, 1, 1);
    t_pointX_t veca, vecb;
    // in [0,1[
    veca.push_back(a1);
    veca.push_back(b1);  // x=t, y=t, z=t
    // in [1,2]
    vecb.push_back(b1);
    vecb.push_back(b1);  // x=(t-1)+1, y=(t-1)+1, z=(t-1)+1
    std::shared_ptr<polynomial_t> pola_ptr =
        std::make_shared<polynomial_t>(veca.begin(), veca.end(), 0, 1);
    std::shared_ptr<polynomial_t> polb_ptr =
        std::make_shared<polynomial_t>(vecb.begin(), vecb.end(), 1, 2);
    piecewise_t pc_C2(pola_ptr);
    pc_C2.add_curve_ptr(polb_ptr);
    // check C0 continuity
    std::string errmsg2(
        "in piecewise polynomial curve test, Error while checking continuity "
        "C0 on ");
    std::string errmsg3(
        "in piecewise polynomial curve test, Error while checking continuity "
        "C1 on ");
    std::string errmsg4(
        "in piecewise polynomial curve test, Error while checking continuity "
        "C2 on ");
    // not C0
    bool isC0 = pc.is_continuous(0);
    if (isC0) {
      std::cout << errmsg2 << " pc " << std::endl;
      error = true;
    }
    // C0
    isC0 = pc_C0.is_continuous(0);
    if (!isC0) {
      std::cout << errmsg2 << " pc_C0 " << std::endl;
      error = true;
    }
    // not C1
    bool isC1 = pc_C0.is_continuous(1);
    if (isC1) {
      std::cout << errmsg3 << " pc_C0 " << std::endl;
      error = true;
    }
    // C1
    isC1 = pc_C1.is_continuous(1);
    if (!isC1) {
      std::cout << errmsg3 << " pc_C1 " << std::endl;
      error = true;
    }
    // not C2
    bool isC2 = pc_C1.is_continuous(2);
    if (isC2) {
      std::cout << errmsg4 << " pc_C1 " << std::endl;
      error = true;
    }
    // C2
    isC2 = pc_C2.is_continuous(2);
    if (!isC2) {
      std::cout << errmsg4 << " pc_C2 " << std::endl;
      error = true;
    }
    // CONVERT PIECEWISE POLYNOMIAL CURVES TO BEZIER AND HERMITE
    std::string errmsg5(
        "in piecewise polynomial curve test, Error while checking piecewise "
        "curve conversion");
    piecewise_t pc_bezier = pc.convert_piecewise_curve_to_bezier<bezier_t>();
    CompareCurves<piecewise_t, piecewise_t>(pc, pc_bezier, errmsg5, error);
    piecewise_t pc_hermite =
        pc.convert_piecewise_curve_to_cubic_hermite<cubic_hermite_spline_t>();
    CompareCurves<piecewise_t, piecewise_t>(pc, pc_hermite, errmsg5, error);
    piecewise_t pc_polynomial_same =
        pc.convert_piecewise_curve_to_polynomial<polynomial_t>();
    CompareCurves<piecewise_t, piecewise_t>(pc, pc_polynomial_same, errmsg5,
                                            error);
    // CONVERT PIECEWISE BEZIER TO POLYNOMIAL AND HERMITE

    std::string errmsg6(
        "in piecewise bezier curve test, Error while checking piecewise curve "
        "conversion");
    piecewise_t pc_bezier1 =
        pc_C0.convert_piecewise_curve_to_bezier<bezier_t>();
    CompareCurves<piecewise_t, piecewise_t>(pc_C0, pc_bezier1, errmsg6, error);
    piecewise_t pc_hermite1 =
        pc_C0
            .convert_piecewise_curve_to_cubic_hermite<cubic_hermite_spline_t>();
    CompareCurves<piecewise_t, piecewise_t>(pc_C0, pc_hermite1, errmsg6, error);
    piecewise_t pc_polynomial1 =
        pc_C0.convert_piecewise_curve_to_polynomial<polynomial_t>();
    CompareCurves<piecewise_t, piecewise_t>(pc_C0, pc_polynomial1, errmsg6,
                                            error);

    // compare compute_derivate and derivate results :

    curve_abc_t* pc_C2_derivate = pc_C2.compute_derivate_ptr(1);
    curve_abc_t* pc_C2_derivate2 = pc_C2.compute_derivate_ptr(2);
    if (pc_C2.min() != pc_C2_derivate->min()) {
      error = true;
      std::cout << "min bounds for curve and it's derivate are not equals."
                << std::endl;
    }
    if (pc_C2.min() != pc_C2_derivate2->min()) {
      error = true;
      std::cout
          << "min bounds for curve and it's second derivate are not equals."
          << std::endl;
    }
    if (pc_C2.max() != pc_C2_derivate->max()) {
      error = true;
      std::cout << "max bounds for curve and it's derivate are not equals."
                << std::endl;
    }
    if (pc_C2.max() != pc_C2_derivate2->max()) {
      error = true;
      std::cout
          << "max bounds for curve and it's second derivate are not equals."
          << std::endl;
    }
    double t = 0.;
    while (t < pc_C2.max()) {
      if (!QuasiEqual(pc_C2.derivate(t, 1), (*pc_C2_derivate)(t))) {
        error = true;
        std::cout << "value not equal between derivate and compute_derivate "
                     "(order 1) at t = "
                  << t << std::endl;
      }
      if (!QuasiEqual(pc_C2.derivate(t, 2), (*pc_C2_derivate2)(t))) {
        error = true;
        std::cout << "value not equal between derivate and compute_derivate "
                     "(order 2) at t = "
                  << t << std::endl;
      }
      t += 0.01;
    }

  } catch (...) {
    error = true;
    std::cout << "Error in piecewiseCurveTest" << std::endl;
  }
}

void curveAbcDimDynamicTest(bool& error) {
  typedef curve_abc<double, double, true> curve_abc_test_t;
  typedef polynomial<double, double, true> polynomial_test_t;
  typedef exact_cubic<double, double, true> exact_cubic_test_t;
  typedef exact_cubic_test_t::spline_constraints spline_constraints_test_t;
  typedef bezier_curve<double, double, true> bezier_test_t;
  typedef cubic_hermite_spline<double, double, true>
      cubic_hermite_spline_test_t;
  curve_abc_test_t* pt_curve_abc;
  // POLYNOMIAL
  point3_t a(1, 1, 1);
  point3_t b(2, 2, 2);
  t_pointX_t vec;
  vec.push_back(a);
  vec.push_back(b);
  polynomial_test_t pol(vec.begin(), vec.end(), 0, 1);
  try {
    pol(0);
    pol(1);
  } catch (...) {
    error = false;
  }
  // BEZIER
  bezier_test_t bc = bezier_from_curve<bezier_test_t>(pol);
  try {
    bc(0);
    bc(1);
  } catch (...) {
    error = false;
  }
  // CUBIC HERMITE
  cubic_hermite_spline_test_t chs =
      hermite_from_curve<cubic_hermite_spline_test_t>(pol);
  try {
    chs(0);
    chs(1);
  } catch (...) {
    error = false;
  }
  // EXACT CUBIC : NOT SUPPORTED, problem to fix later
  ndcurves::T_Waypoint waypoints;
  for (double i = 0; i <= 1; i = i + 0.2) {
    waypoints.push_back(std::make_pair(i, point3_t(i, i, i)));
  }
  std::string errmsg(
      "Error in ExactCubicVelocityConstraintsTest (1); while checking that "
      "given wayPoints are crossed (expected / "
      "obtained)");
  spline_constraints_test_t constraints(3);
  constraints.end_vel = point3_t(0, 0, 0);
  constraints.init_vel = point3_t(0, 0, 0);
  constraints.end_acc = point3_t(0, 0, 0);
  constraints.init_acc = point3_t(0, 0, 0);
  exact_cubic_test_t ec(waypoints.begin(), waypoints.end(), constraints);
  try {
    ec(0);
    ec(1);
  } catch (...) {
    error = false;
  }
  // Test with pointer to curve_abc type
  try {
    pt_curve_abc = &pol;
    (*pt_curve_abc)(0);
    (*pt_curve_abc)(1);
    pt_curve_abc = &bc;
    (*pt_curve_abc)(0);
    (*pt_curve_abc)(1);
    pt_curve_abc = &chs;
    (*pt_curve_abc)(0);
    (*pt_curve_abc)(1);
    pt_curve_abc = &ec;
    (*pt_curve_abc)(0);
    (*pt_curve_abc)(1);
  } catch (...) {
    error = false;
  }
}

void PiecewisePolynomialCurveFromDiscretePoints(bool& error) {
  std::string errMsg(
      "PiecewisePolynomialCurveFromDiscretePoints, Error, value on curve is "
      "wrong : ");
  point3_t p0(0., 0., 0.);
  point3_t p1(1., 2., 3.);
  point3_t p2(4., 4., 4.);
  point3_t p3(10., 10., 10.);
  point3_t d0(1., 1., 1.);
  point3_t d1(2., 2., 2.);
  point3_t d2(3., 3., 3.);
  point3_t d3(5., 5., 5.);
  point3_t dd0(1.5, 1.5, 1.5);
  point3_t dd1(2.5, 2.5, 2.5);
  point3_t dd2(3.5, 3.5, 3.5);
  point3_t dd3(5.5, 5.5, 5.5);
  double t0 = 1.0;
  double t1 = 1.5;
  double t2 = 3.0;
  double t3 = 10.0;
  t_pointX_t points;
  points.push_back(p0);
  points.push_back(p1);
  points.push_back(p2);
  points.push_back(p3);
  t_pointX_t points_derivative;
  points_derivative.push_back(d0);
  points_derivative.push_back(d1);
  points_derivative.push_back(d2);
  points_derivative.push_back(d3);
  t_pointX_t points_second_derivative;
  points_second_derivative.push_back(dd0);
  points_second_derivative.push_back(dd1);
  points_second_derivative.push_back(dd2);
  points_second_derivative.push_back(dd3);
  std::vector<double> time_points;
  time_points.push_back(t0);
  time_points.push_back(t1);
  time_points.push_back(t2);
  time_points.push_back(t3);

  // Piecewise polynomial curve C0 => Linear interpolation between points
  piecewise_t ppc_C0 =
      piecewise_t::convert_discrete_points_to_polynomial<polynomial_t>(
          points, time_points);
  if (!ppc_C0.is_continuous(0)) {
    std::cout << "PiecewisePolynomialCurveFromDiscretePoints, Error, piecewise "
                 "curve is not C0"
              << std::endl;
    error = true;
  }
  for (std::size_t i = 0; i < points.size(); i++) {
    ComparePoints(points[i], ppc_C0(time_points[i]), errMsg, error);
  }
  point3_t pos_between_po_and_p1((p1[0] + p0[0]) / 2.0, (p1[1] + p0[1]) / 2.0,
                                 (p1[2] + p0[2]) / 2.0);
  double time_between_po_and_p1 = (t0 + t1) / 2.0;
  ComparePoints(pos_between_po_and_p1, ppc_C0(time_between_po_and_p1), errMsg,
                error);

  // Piecewise polynomial curve C1
  piecewise_t ppc_C1 =
      piecewise_t::convert_discrete_points_to_polynomial<polynomial_t>(
          points, points_derivative, time_points);
  if (!ppc_C1.is_continuous(1)) {
    std::cout << "PiecewisePolynomialCurveFromDiscretePoints, Error, piecewise "
                 "curve is not C1"
              << std::endl;
    error = true;
  }
  for (std::size_t i = 0; i < points.size(); i++) {
    ComparePoints(points[i], ppc_C1(time_points[i]), errMsg, error);
    ComparePoints(points_derivative[i], ppc_C1.derivate(time_points[i], 1),
                  errMsg, error);
  }

  // Piecewise polynomial curve C2
  piecewise_t ppc_C2 =
      piecewise_t::convert_discrete_points_to_polynomial<polynomial_t>(
          points, points_derivative, points_second_derivative, time_points);
  if (!ppc_C2.is_continuous(2)) {
    std::cout << "PiecewisePolynomialCurveFromDiscretePoints, Error, piecewise "
                 "curve is not C1"
              << std::endl;
    error = true;
  }
  for (std::size_t i = 0; i < points.size(); i++) {
    ComparePoints(points[i], ppc_C2(time_points[i]), errMsg, error);
    ComparePoints(points_derivative[i], ppc_C2.derivate(time_points[i], 1),
                  errMsg, error);
    ComparePoints(points_second_derivative[i],
                  ppc_C2.derivate(time_points[i], 2), errMsg, error);
  }
}

void PiecewisePolynomialCurveFromFile(bool& error) {
  std::string filename_pos(TEST_DATA_PATH "discrete_points_pos.txt");
  std::string filename_vel(TEST_DATA_PATH "discrete_points_vel.txt");
  std::string filename_acc(TEST_DATA_PATH "discrete_points_acc.txt");
  std::string filename_error(TEST_DATA_PATH "discrete_points_error.txt");

  piecewise_t c_pos = piecewise_t::load_piecewise_from_text_file<polynomial_t>(
      filename_pos, 0.01, 3);
  if (c_pos.min() != 0.) {
    std::cout << "PiecewisePolynomialCurveFromFile, Error, t_min should be 0"
              << std::endl;
    error = true;
  }
  if (c_pos.max() != 0.03) {
    std::cout << "PiecewisePolynomialCurveFromFile, Error, t_max should be 0.03"
              << std::endl;
    error = true;
  }
  pointX_t p0(3), p2(3);
  p0 << -0.003860389372941039, 0.0012353625242474164, 0.009005041639999767;
  p2 << -0.0028803627898293283, 0.0011918668401150736, 0.009005041639999767;
  if (!c_pos(0.).isApprox(p0)) {
    std::cout << "PiecewisePolynomialCurveFromFile, Error, points do not match"
              << std::endl;
    error = true;
  }
  if (!c_pos(0.02).isApprox(p2)) {
    std::cout << "PiecewisePolynomialCurveFromFile, Error, points do not match"
              << std::endl;
    error = true;
  }

  piecewise_t c_vel = piecewise_t::load_piecewise_from_text_file<polynomial_t>(
      filename_vel, 0.05, 3);
  if (c_pos.min() != 0.) {
    std::cout << "PiecewisePolynomialCurveFromFile, Error, t_min should be 0"
              << std::endl;
    error = true;
  }
  if (!QuasiEqual(c_vel.max(), 0.15)) {
    std::cout << "PiecewisePolynomialCurveFromFile, Error, t_max should be 0.15"
              << std::endl;
    error = true;
  }
  pointX_t p3(3);
  p3 << 0.2968141884672718, 0.0012916907964522569, 0.00951023474821927;
  if (!c_vel(0.).isApprox(p0)) {
    std::cout << "PiecewisePolynomialCurveFromFile, Error, points do not match"
              << std::endl;
    error = true;
  }
  if (!c_vel(0.15).isApprox(p3)) {
    std::cout << "PiecewisePolynomialCurveFromFile, Error, points do not match"
              << std::endl;
    error = true;
  }
  if (!c_vel.derivate(0., 1).isZero()) {
    std::cout << "PiecewisePolynomialCurveFromFile, Error, c_vel derivative at "
                 "0. should be null"
              << std::endl;
    error = true;
  }
  if (!c_vel.derivate(0.1, 1).isZero()) {
    std::cout << "PiecewisePolynomialCurveFromFile, Error, c_vel derivative at "
                 "0.1 should be null"
              << std::endl;
    error = true;
  }

  piecewise_t c_acc = piecewise_t::load_piecewise_from_text_file<polynomial_t>(
      filename_acc, 0.001, 3);
  if (c_acc.min() != 0.) {
    std::cout << "PiecewisePolynomialCurveFromFile, Error, t_min should be 0"
              << std::endl;
    error = true;
  }
  if (!QuasiEqual(c_acc.max(), 7.85)) {
    std::cout << "PiecewisePolynomialCurveFromFile, Error, t_max should be 7.85"
              << std::endl;
    error = true;
  }
  if (!c_acc(0.).isApprox(p0)) {
    std::cout << "PiecewisePolynomialCurveFromFile, Error, points do not match"
              << std::endl;
    error = true;
  }
  pointX_t p5200(3);
  p5200 << 0.30273356072723845, -0.07619420199174821, 0.010015348526727433;
  if (!c_acc(5.2).isApprox(p5200)) {
    std::cout << "PiecewisePolynomialCurveFromFile, Error, points do not match"
              << std::endl;
    error = true;
  }
  if (!c_acc.derivate(0., 1).isZero()) {
    std::cout << "PiecewisePolynomialCurveFromFile, Error, c_acc derivative at "
                 "0 should be null"
              << std::endl;
    error = true;
  }
  if (!c_acc.derivate(0.5, 1).isZero()) {
    std::cout << "PiecewisePolynomialCurveFromFile, Error, c_acc derivative "
                 "should at 0.5 be null"
              << std::endl;
    error = true;
  }
  if (!c_acc.derivate(0., 2).isZero()) {
    std::cout << "PiecewisePolynomialCurveFromFile, Error, c_acc second "
                 "derivative at 0 should be null"
              << std::endl;
    error = true;
  }
  if (!c_acc.derivate(5., 2).isZero()) {
    std::cout << "PiecewisePolynomialCurveFromFile, Error, c_acc second "
                 "derivative at 5 should be null"
              << std::endl;
    error = true;
  }

  try {
    piecewise_t c_error =
        piecewise_t::load_piecewise_from_text_file<polynomial_t>(filename_acc,
                                                                 0.01, 4);
    std::cout << "PiecewisePolynomialCurveFromFile, Error, dimension do not "
                 "match, an error should be raised"
              << std::endl;
    error = true;
  } catch (std::invalid_argument& /*e*/) {
  }
  try {
    piecewise_t c_error =
        piecewise_t::load_piecewise_from_text_file<polynomial_t>(filename_error,
                                                                 0.01, 3);
    std::cout << "PiecewisePolynomialCurveFromFile, Error, "
                 "discrete_points_error should not be parsed correctly"
              << std::endl;
    error = true;
  } catch (std::invalid_argument& /*e*/) {
  }
}

void serializationCurvesTest(bool& error) {
  try {
    std::string errMsg1(
        "in serializationCurveTest, Error While serializing Polynomial : ");
    std::string errMsg2(
        "in serializationCurveTest, Error While serializing Bezier : ");
    std::string errMsg3(
        "in serializationCurveTest, Error While serializing Cubic Hermite : ");
    std::string errMsg4(
        "in serializationCurveTest, Error While serializing Piecewise curves "
        ": ");
    std::string errMsg5(
        "in serializationCurveTest, Error While serializing Exact cubic : ");
    std::string errMsg6(
        "in serializationCurveTest, Error While serializing using abstract "
        "pointers : ");
    point3_t a(1, 1, 1);  // in [0,1[
    point3_t b(2, 1, 1);  // in [1,2[
    point3_t c(3, 1, 1);  // in [2,3]
    point3_t res;
    t_pointX_t vec1, vec2, vec3;
    vec1.push_back(a);  // x=1, y=1, z=1
    vec2.push_back(b);  // x=2, y=1, z=1
    vec3.push_back(c);  // x=3, y=1, z=1
    polynomial_t pol1(vec1.begin(), vec1.end(), 0, 1);
    polynomial_t pol2(vec2.begin(), vec2.end(), 1, 2);
    polynomial_t pol3(vec3.begin(), vec3.end(), 2, 3);
    piecewise_t ppc;
    ppc.add_curve<polynomial_t>(pol1);
    ppc.add_curve<polynomial_t>(pol2);
    ppc.add_curve<polynomial_t>(pol3);
    std::string fileName("fileTest.test");
    std::string fileName1("fileTest1.test");
    // Simple curves
    // Test serialization on Polynomial
    pol1.saveAsText<polynomial_t>(fileName1);
    polynomial_t pol_test;
    pol_test.loadFromText<polynomial_t>(fileName1);
    CompareCurves<polynomial_t, polynomial_t>(pol1, pol_test, errMsg1, error);
    // Test serialization on Bezier
    bezier_t bc = bezier_from_curve<bezier_t>(pol1);
    bc.saveAsText<bezier_t>(fileName);
    bezier_t bc_test;
    bc_test.loadFromText<bezier_t>(fileName);
    CompareCurves<polynomial_t, bezier_t>(pol1, bc_test, errMsg2, error);
    // Test serialization on Cubic Hermite
    cubic_hermite_spline_t chs =
        hermite_from_curve<cubic_hermite_spline_t>(pol1);
    chs.saveAsText<cubic_hermite_spline_t>(fileName);
    cubic_hermite_spline_t chs_test;
    chs_test.loadFromText<cubic_hermite_spline_t>(fileName);
    CompareCurves<polynomial_t, cubic_hermite_spline_t>(pol1, chs_test, errMsg3,
                                                        error);
    // Piecewise curves
    // Test serialization on Piecewise Polynomial curve
    ppc.saveAsText<piecewise_t>(fileName);
    piecewise_t ppc_test, ppc_test_binary;
    ppc_test.loadFromText<piecewise_t>(fileName);
    CompareCurves<piecewise_t, piecewise_t>(ppc, ppc_test, errMsg4, error);
    ppc.saveAsBinary<piecewise_t>(fileName);
    ppc_test_binary.loadFromBinary<piecewise_t>(fileName);
    CompareCurves<piecewise_t, piecewise_t>(ppc, ppc_test_binary, errMsg4,
                                            error);

    // Test serialization on Piecewise Bezier curve
    piecewise_t pbc = ppc.convert_piecewise_curve_to_bezier<bezier_t>();
    pbc.saveAsText<piecewise_t>(fileName);
    piecewise_t pbc_test;
    pbc_test.loadFromText<piecewise_t>(fileName);
    CompareCurves<piecewise_t, piecewise_t>(ppc, pbc_test, errMsg4, error);
    // Test serialization on Piecewise Cubic Hermite curve
    piecewise_t pchc =
        ppc.convert_piecewise_curve_to_cubic_hermite<cubic_hermite_spline_t>();
    pchc.saveAsText<piecewise_t>(fileName);
    piecewise_t pchc_test;
    pchc_test.loadFromText<piecewise_t>(fileName);
    CompareCurves<piecewise_t, piecewise_t>(ppc, pchc_test, errMsg4, error);
    // Test serialization on exact cubic
    ndcurves::T_Waypoint waypoints;
    for (double i = 0; i <= 1; i = i + 0.2) {
      waypoints.push_back(std::make_pair(i, point3_t(i, i, i)));
    }
    spline_constraints_t constraints(3);
    constraints.end_vel = point3_t(0.1, 0, 0);
    constraints.init_vel = point3_t(0.2, 0, 0);
    constraints.end_acc = point3_t(0.01, 0, 0);
    constraints.init_acc = point3_t(0.01, 0, 0);
    exact_cubic_t ec(waypoints.begin(), waypoints.end(), constraints);
    ec.saveAsText<exact_cubic_t>(fileName);
    exact_cubic_t ec_test;
    ec_test.loadFromText<exact_cubic_t>(fileName);
    CompareCurves<exact_cubic_t, exact_cubic_t>(ec, ec_test, errMsg5, error);
    // Test with pointer on abstract struct curve_abc
    // Polynomial
    curve_abc_t* pt_0;
    curve_abc_t* pt_1;
    pol_test = polynomial_t();
    pt_0 = &pol1;
    pt_1 = &pol_test;
    (*pt_0).saveAsText<polynomial_t>(fileName);
    (*pt_1).loadFromText<polynomial_t>(fileName);
    CompareCurves<polynomial_t, polynomial_t>(
        pol1, (*dynamic_cast<polynomial_t*>(pt_1)), errMsg6, error);
    // Piecewise Polynomial
    pt_0 = NULL;
    pt_1 = NULL;
    ppc_test = piecewise_t();
    pt_0 = &ppc;
    pt_1 = &ppc_test;
    (*pt_0).saveAsText<piecewise_t>(fileName);
    (*pt_1).loadFromText<piecewise_t>(fileName);
    CompareCurves<piecewise_t, piecewise_t>(
        ppc, (*dynamic_cast<piecewise_t*>(pt_1)), errMsg6, error);
  } catch (...) {
    error = true;
    std::cout << "Error in serializationCurvesTest" << std::endl;
  }
}

void polynomialFromBoundaryConditions(bool& error) {
  pointX_t zeros = point3_t(0., 0., 0.);
  pointX_t p0 = point3_t(0., 1., 0.);
  pointX_t p1 = point3_t(1., 2., -3.);
  pointX_t dp0 = point3_t(-8., 4., 6.);
  pointX_t dp1 = point3_t(10., -10., 10.);
  pointX_t ddp0 = point3_t(-1., 7., 4.);
  pointX_t ddp1 = point3_t(12., -8., 2.5);
  double min = 0.5;
  double max = 2.;
  // C0 : order 1
  polynomial_t polC0 = polynomial_t(p0, p1, min, max);
  if (polC0.min() != min) {
    error = true;
    std::cout
        << "polynomialFromBoundaryConditions C0: min interval not respected."
        << std::endl;
  }
  if (polC0.max() != max) {
    error = true;
    std::cout
        << "polynomialFromBoundaryConditions C0: max interval not respected."
        << std::endl;
  }
  if (polC0(min) != p0) {
    error = true;
    std::cout
        << "polynomialFromBoundaryConditions C0: initial value not respected"
        << std::endl;
  }
  if (polC0(max) != p1) {
    error = true;
    std::cout
        << "polynomialFromBoundaryConditions C0: final value not respected"
        << std::endl;
  }
  if (polC0.degree_ != 1) {
    error = true;
    std::cout << "polynomialFromBoundaryConditions C0: curve is not degree 1 "
              << std::endl;
  }
  if (polC0((max + min) / 2.) != (p0 * 0.5 + p1 * 0.5)) {
    error = true;
    std::cout << "polynomialFromBoundaryConditions C0: middle point doesn't "
                 "have the right value' "
              << std::endl;
  }
  // C1 : order 3
  polynomial_t polC1 = polynomial_t(p0, dp0, p1, dp1, min, max);
  if (polC1.min() != min) {
    error = true;
    std::cout
        << "polynomialFromBoundaryConditions C1: min interval not respected."
        << std::endl;
  }
  if (polC1.max() != max) {
    error = true;
    std::cout
        << "polynomialFromBoundaryConditions C1: max interval not respected."
        << std::endl;
  }
  if (polC1(min) != p0) {
    error = true;
    std::cout
        << "polynomialFromBoundaryConditions C1: initial value not respected"
        << std::endl;
  }
  if (!QuasiEqual(polC1(max), p1)) {
    error = true;
    std::cout
        << "polynomialFromBoundaryConditions C1: final value not respected"
        << std::endl;
    std::cout << "p1 = " << p1.transpose()
              << " curve end = " << polC1(max).transpose() << std::endl;
  }
  if (polC1.derivate(min, 1) != dp0) {
    error = true;
    std::cout << "polynomialFromBoundaryConditions C1: initial derivative "
                 "value not respected"
              << std::endl;
  }
  if (!QuasiEqual(polC1.derivate(max, 1), dp1)) {
    error = true;
    std::cout << "polynomialFromBoundaryConditions C1: final derivative value "
                 "not respected"
              << std::endl;
    std::cout << "dp1 = " << dp1.transpose() << " curve end derivative = "
              << polC1.derivate(max, 1).transpose() << std::endl;
  }
  if (polC1.degree_ != 3) {
    error = true;
    std::cout << "polynomialFromBoundaryConditions C1: curve is not degree 3 "
              << std::endl;
  }
  // C2 : order 5
  polynomial_t polC2 = polynomial_t(p0, dp0, ddp0, p1, dp1, ddp1, min, max);
  if (polC2.min() != min) {
    error = true;
    std::cout
        << "polynomialFromBoundaryConditions C2: min interval not respected."
        << std::endl;
  }
  if (polC2.max() != max) {
    error = true;
    std::cout
        << "polynomialFromBoundaryConditions C2: max interval not respected."
        << std::endl;
  }
  if (polC2(min) != p0) {
    error = true;
    std::cout
        << "polynomialFromBoundaryConditions C2: initial value not respected"
        << std::endl;
  }
  if (!QuasiEqual(polC2(max), p1)) {
    error = true;
    std::cout
        << "polynomialFromBoundaryConditions C2: final value not respected"
        << std::endl;
  }
  if (polC2.derivate(min, 1) != dp0) {
    error = true;
    std::cout << "polynomialFromBoundaryConditions C2: initial derivative "
                 "value not respected"
              << std::endl;
  }
  if (!QuasiEqual(polC2.derivate(max, 1), dp1)) {
    error = true;
    std::cout << "polynomialFromBoundaryConditions C2: final derivative value "
                 "not respected"
              << std::endl;
  }
  if (polC2.derivate(min, 2) != ddp0) {
    error = true;
    std::cout << "polynomialFromBoundaryConditions C2: initial second "
                 "derivative value not respected"
              << std::endl;
  }
  if (!QuasiEqual(polC2.derivate(max, 2), ddp1)) {
    error = true;
    std::cout << "polynomialFromBoundaryConditions C2: final second derivative "
                 "value not respected"
              << std::endl;
  }
  if (polC2.degree_ != 5) {
    error = true;
    std::cout << "polynomialFromBoundaryConditions C2: curve is not degree 5 "
              << std::endl;
  }
  // check if the exeptions are correctly raised :
  try {
    polynomial_t polC0Err = polynomial_t(p0, p1, max, min);
    error = true;
    std::cout << "Created a polynomial with tMin > tMax without error. "
              << std::endl;
  } catch (const invalid_argument& /*e*/) {
  }
  try {
    polynomial_t polC1Err = polynomial_t(p0, dp0, p1, dp1, max, min);
    error = true;
    std::cout << "Created a polynomial with tMin > tMax without error. "
              << std::endl;
  } catch (const invalid_argument& /*e*/) {
  }
  try {
    polynomial_t polC2Err =
        polynomial_t(p0, dp0, ddp0, p1, dp1, ddp1, max, min);
    error = true;
    std::cout << "Created a polynomial with tMin > tMax without error. "
              << std::endl;
  } catch (const invalid_argument& /*e*/) {
  }
}

void so3LinearTest(bool& error) {
  quaternion_t q0(1, 0, 0, 0);
  quaternion_t q1(0.7071, 0.7071, 0, 0);
  const double tMin = 0.;
  const double tMax = 1.5;
  SO3Linear_t so3Traj(q0, q1, tMin, tMax);

  if (so3Traj.min() != tMin) {
    error = true;
    std::cout << "Min bound not respected" << std::endl;
  }
  if (so3Traj.max() != tMax) {
    error = true;
    std::cout << "Max bound not respected" << std::endl;
  }
  if (!so3Traj.computeAsQuaternion(tMin).isApprox(q0)) {
    error = true;
    std::cout << "evaluate at t=0 is not the init quaternion" << std::endl;
  }
  if (so3Traj(tMin) != q0.toRotationMatrix()) {
    error = true;
    std::cout << "evaluate at t=0 is not the init rotation" << std::endl;
  }
  if (!so3Traj.computeAsQuaternion(tMax).isApprox(q1)) {
    error = true;
    std::cout << "evaluate at t=max is not the final quaternion" << std::endl;
  }
  if (so3Traj(tMax) != q1.toRotationMatrix()) {
    error = true;
    std::cout << "evaluate at t=max is not the final rotation" << std::endl;
  }
  // check derivatives :
  if (so3Traj.derivate(tMin, 1) != so3Traj.derivate(1., 1)) {
    error = true;
    std::cout << "first order derivative should be constant." << std::endl;
  }
  if (so3Traj.derivate(tMin, 2) != point3_t::Zero(3)) {
    error = true;
    std::cout << "second order derivative should be null" << std::endl;
  }
  point3_t angular_vel = so3Traj.derivate(tMin, 1);
  if (angular_vel[1] != 0. || angular_vel[2] != 0) {
    error = true;
    std::cout << "Angular velocity around y and z axis should be null"
              << std::endl;
  }

  constant3_t so3Derivate1 = so3Traj.compute_derivate(1);
  if (so3Derivate1(1.) != so3Traj.derivate(1., 1)) {
    error = true;
    std::cout << "compute_derivate curve do not equal derivate call"
              << std::endl;
  }

  constant3_t so3Derivate2 = so3Traj.compute_derivate(2);
  if (so3Derivate2(1.) != point3_t::Zero(3)) {
    error = true;
    std::cout << "compute_derivate curve do not equal derivate call"
              << std::endl;
  }

  // check if errors are correctly raised :
  try {
    so3Traj(-0.1);
    error = true;
    std::cout << "SO3Linear: calling () with t < tmin should raise an "
                 "invalid_argument error"
              << std::endl;
  } catch (std::invalid_argument& /*e*/) {
  }
  try {
    so3Traj(1.7);
    error = true;
    std::cout << "SO3Linear: calling () with t > tmin should raise an "
                 "invalid_argument error"
              << std::endl;
  } catch (std::invalid_argument& /*e*/) {
  }
  try {
    so3Traj.derivate(0, 0);
    error = true;
    std::cout << "SO3Linear: calling derivate with order = 0 should raise an "
                 "invalid_argument error"
              << std::endl;
  } catch (std::invalid_argument& /*e*/) {
  }

  SO3Linear_t so3TrajMatrix(q0.toRotationMatrix(), q1.toRotationMatrix(), tMin,
                            tMax);
  std::string errmsg(
      "SO3Linear built from quaternion or from matrix are not identical.");
  CompareCurves(so3Traj, so3TrajMatrix, errmsg, error, 1e-3);
}

void SO3serializationTest(bool& error) {
  std::string fileName("fileTest");
  std::string errmsg(
      "SO3serializationTest : curve serialized is not equivalent to the "
      "original curve.");
  quaternion_t q0(0.544, -0.002, -0.796, 0.265);
  quaternion_t q1(0.7071, 0.7071, 0, 0);
  q1.normalize();
  q0.normalize();
  SO3Linear_t so3Traj(q0, q1, 0.5, 2.2);

  so3Traj.saveAsText<SO3Linear_t>(fileName + ".txt");
  SO3Linear_t so3_from_txt;
  so3_from_txt.loadFromText<SO3Linear_t>(fileName + ".txt");
  CompareCurves<SO3Linear_t, SO3Linear_t>(
      so3Traj, so3_from_txt, errmsg + " For text serialization", error);

  so3Traj.saveAsXML<SO3Linear_t>(fileName + ".xml", "so3Curve");
  SO3Linear_t so3_from_xml;
  so3_from_xml.loadFromXML<SO3Linear_t>(fileName + ".xml", "so3Curve");
  CompareCurves<SO3Linear_t, SO3Linear_t>(
      so3Traj, so3_from_xml, errmsg + " For XML serialization", error);

  so3Traj.saveAsBinary<SO3Linear_t>(fileName);
  SO3Linear_t so3_from_binary;
  so3_from_binary.loadFromBinary<SO3Linear_t>(fileName);
  CompareCurves<SO3Linear_t, SO3Linear_t>(
      so3Traj, so3_from_binary, errmsg + " For binary serialization", error);
}

void se3CurveTest(bool& error) {
  quaternion_t q0(1, 0, 0, 0);
  quaternion_t q1(0., 1., 0, 0);
  pointX_t p0 = point3_t(1., 1.5, -2.);
  pointX_t p1 = point3_t(3., 0, 1.);

  double min = 0.5, max = 2.;

  // constructor from init/end position/rotation : automatically create a linear
  // interpolation for position and slerp for rotation
  SE3Curve_t cLinear(p0, p1, q0, q1, min, max);
  transform_t transformInit = cLinear(min);
  transform_t transformEnd = cLinear(max);
  transform_t transformMid = cLinear((max + min) / 2.);
  if (!transformInit.translation().isApprox(p0)) {
    error = true;
    std::cout << "Init position of the curve is not correct." << std::endl;
  }
  if (!transformInit.rotation().isApprox(q0.toRotationMatrix())) {
    error = true;
    std::cout << "Init rotation of the curve is not correct." << std::endl;
  }
  if (!transformEnd.translation().isApprox(p1)) {
    error = true;
    std::cout << "End position of the curve is not correct." << std::endl;
  }
  if (!transformEnd.rotation().isApprox(q1.toRotationMatrix())) {
    error = true;
    std::cout << "End rotation of the curve is not correct." << std::endl;
  }
  quaternion_t qMid(sqrt(2.) / 2., sqrt(2.) / 2., 0, 0);
  point3_t pMid = (p0 + p1) / 2.;
  if (!transformMid.translation().isApprox(pMid)) {
    error = true;
    std::cout << "Mid position of the curve is not correct." << std::endl;
  }
  if (!transformMid.rotation().isApprox(qMid.toRotationMatrix())) {
    error = true;
    std::cout << "Mid rotation of the curve is not correct." << std::endl;
  }

  // constructor with specific translation curve
  SE3Curve_t cBezier;
  {  // inner scope to check what happen when translation_bezier is out of scope
    point3_t a(1, 2, 3);
    point3_t b(2, 3, 4);
    point3_t c(3, 4, 5);
    point3_t d(3, 6, 7);
    std::vector<point3_t> params;
    params.push_back(a);
    params.push_back(b);
    params.push_back(c);
    params.push_back(d);
    std::shared_ptr<bezier3_t> translation_bezier =
        std::make_shared<bezier3_t>(params.begin(), params.end(), min, max);
    cBezier = SE3Curve_t(translation_bezier, q0.toRotationMatrix(),
                         q1.toRotationMatrix());
    p0 = (*translation_bezier)(min);
    p1 = (*translation_bezier)(max);
    pMid = (*translation_bezier)((max + min) / 2.);
    if (cBezier.min() != min) {
      error = true;
      std::cout << "SE3 constructor from translation bezier do not respect the "
                   "min time interval"
                << std::endl;
    }
    if (cBezier.max() != max) {
      error = true;
      std::cout << "SE3 constructor from translation bezier do not respect the "
                   "max time interval"
                << std::endl;
    }
    double t = min;
    while (t < max) {
      if (!cBezier(t).translation().isApprox((*translation_bezier)(t))) {
        error = true;
        std::cout << "SE3 translation is not equivalent to bezier for t = " << t
                  << std::endl;
      }
      t += 0.1;
    }
    // check the derivatives for translation:
    for (size_t i = 1; i < 3; i++) {
      t = min;
      while (t < max) {
        if (!cBezier.derivate(t, i).head<3>().isApprox(
                translation_bezier->derivate(t, i))) {
          error = true;
          std::cout
              << "SE3 curve derivative is not equivalent to bezier for t = "
              << t << " and order = " << i << std::endl;
        }
        t += 0.1;
      }
    }
  }
  // check the rotation
  transformInit = cBezier(min);
  transformEnd = cBezier(max);
  transformMid = cBezier((max + min) / 2.);
  if (!transformInit.translation().isApprox(p0)) {
    error = true;
    std::cout << "Init position of the curve is not correct." << std::endl;
  }
  if (!transformInit.rotation().isApprox(q0.toRotationMatrix())) {
    error = true;
    std::cout << "Init rotation of the curve is not correct." << std::endl;
  }
  if (!transformEnd.translation().isApprox(p1)) {
    error = true;
    std::cout << "End position of the curve is not correct." << std::endl;
  }
  if (!transformEnd.rotation().isApprox(q1.toRotationMatrix())) {
    error = true;
    std::cout << "End rotation of the curve is not correct." << std::endl;
  }
  if (!transformMid.translation().isApprox(pMid)) {
    error = true;
    std::cout << "Mid position of the curve is not correct." << std::endl;
  }
  if (!transformMid.rotation().isApprox(qMid.toRotationMatrix())) {
    error = true;
    std::cout << "Mid rotation of the curve is not correct." << std::endl;
  }

  // check derivatives for rotation:
  if (cBezier.derivate(min, 1).tail<3>() !=
      cBezier.derivate(max, 1).tail<3>()) {
    error = true;
    std::cout
        << "SE3 curve : first order derivative for rotation should be constant."
        << std::endl;
  }
  if (cBezier.derivate(min, 2).tail<3>() != point3_t::Zero(3)) {
    error = true;
    std::cout
        << "SE3 curve : second order derivative for rotation should be null"
        << std::endl;
  }

  // check accessor to translation curves :
  curve_translation_ptr_t translation = cBezier.translation_curve();
  if (translation->operator()(min) != cBezier(min).translation()) {
    error = true;
    std::cout << "SE3 curve : translation curve not equal to se3.translation"
              << std::endl;
  }
  if (translation->operator()(max) != cBezier(max).translation()) {
    error = true;
    std::cout << "SE3 curve : translation curve not equal to se3.translation"
              << std::endl;
  }
  if (translation->operator()((max + min) / 2.) !=
      cBezier((max + min) / 2.).translation()) {
    error = true;
    std::cout << "SE3 curve : translation curve not equal to se3.translation"
              << std::endl;
  }
  // check accessor to rotation curves :
  curve_rotation_ptr_t rotation = cBezier.rotation_curve();
  if (!rotation->operator()(min).isApprox(cBezier(min).rotation())) {
    error = true;
    std::cout << "SE3 curve : rotation curve not equal to se3.rotation"
              << std::endl;
  }
  if (!rotation->operator()(max).isApprox(cBezier(max).rotation())) {
    error = true;
    std::cout << "SE3 curve : rotation curve not equal to se3.rotation"
              << std::endl;
  }
  if (!rotation->operator()((max + min) / 2.)
           .isApprox(cBezier((max + min) / 2.).rotation())) {
    error = true;
    std::cout << "SE3 curve : rotation curve not equal to se3.rotation"
              << std::endl;
  }

  // check if errors are correctly raised
  try {
    cBezier(0.1);
    error = true;
    std::cout << "SE3 curve: calling () with t < tmin should raise an "
                 "invalid_argument error"
              << std::endl;
  } catch (std::invalid_argument& /*e*/) {
  }
  try {
    cBezier(2.3);
    error = true;
    std::cout << "SE3 curve: calling () with t > tmin should raise an "
                 "invalid_argument error"
              << std::endl;
  } catch (std::invalid_argument& /*e*/) {
  }
  try {
    cBezier.derivate(0.6, 0);
    error = true;
    std::cout << "SE3 curve: calling derivate with order = 0 should raise an "
                 "invalid_argument error"
              << std::endl;
  } catch (std::invalid_argument& /*e*/) {
  }
}

void Se3serializationTest(bool& error) {
  std::string fileName("fileTest");
  std::string errmsg(
      "SE3serializationTest : curve serialized is not equivalent to the "
      "original curve.");
  quaternion_t q0(1, 0, 0, 0);
  quaternion_t q1(0., 1., 0, 0);
  pointX_t p0 = point3_t(1., 1.5, -2.);
  pointX_t p1 = point3_t(3., 0, 1.);
  double min = 0.5, max = 2.;
  // constructor from init/end position/rotation : automatically create a linear
  // interpolation for position and slerp for rotation
  SE3Curve_t cLinear(p0, p1, q0, q1, min, max);

  cLinear.saveAsText<SE3Curve_t>(fileName + ".txt");
  SE3Curve_t se3_from_txt;
  se3_from_txt.loadFromText<SE3Curve_t>(fileName + ".txt");
  CompareCurves<SE3Curve_t, SE3Curve_t>(
      cLinear, se3_from_txt, errmsg + " For text serialization", error);

  cLinear.saveAsXML<SE3Curve_t>(fileName + ".xml", "se3Curve");
  SE3Curve_t se3_from_xml;
  se3_from_xml.loadFromXML<SE3Curve_t>(fileName + ".xml", "se3Curve");
  CompareCurves<SE3Curve_t, SE3Curve_t>(
      cLinear, se3_from_xml, errmsg + " For XML serialization", error);

  cLinear.saveAsBinary<SE3Curve_t>(fileName);
  SE3Curve_t se3_from_binary;
  se3_from_binary.loadFromBinary<SE3Curve_t>(fileName);
  CompareCurves<SE3Curve_t, SE3Curve_t>(
      cLinear, se3_from_binary, errmsg + " For binary serialization", error);

  // constructor with specific translation curve
  SE3Curve_t cBezier;
  {  // inner scope to check what happen when translation_bezier is out of scope
    quaternion_t q0(0.544, -0.002, -0.796, 0.265);
    quaternion_t q1(0.7071, 0.7071, 0, 0);
    q1.normalize();
    q0.normalize();
    point3_t a(1, 2, 3);
    point3_t b(2, 3, 4);
    point3_t c(3, 4, 5);
    point3_t d(3, 6, 7);
    std::vector<point3_t> params;
    params.push_back(a);
    params.push_back(b);
    params.push_back(c);
    params.push_back(d);
    std::shared_ptr<bezier3_t> translation_bezier =
        std::make_shared<bezier3_t>(params.begin(), params.end(), min, max);
    cBezier = SE3Curve_t(translation_bezier, q0, q1);
  }

  cBezier.saveAsText<SE3Curve_t>(fileName + ".txt");
  SE3Curve_t se3_from_txt_bezier;
  se3_from_txt_bezier.loadFromText<SE3Curve_t>(fileName + ".txt");
  CompareCurves<SE3Curve_t, SE3Curve_t>(
      cBezier, se3_from_txt_bezier, errmsg + " For text serialization", error);

  cBezier.saveAsXML<SE3Curve_t>(fileName + ".xml", "se3Curve");
  SE3Curve_t se3_from_xml_bezier;
  se3_from_xml_bezier.loadFromXML<SE3Curve_t>(fileName + ".xml", "se3Curve");
  CompareCurves<SE3Curve_t, SE3Curve_t>(
      cBezier, se3_from_xml_bezier, errmsg + " For XML serialization", error);

  cBezier.saveAsBinary<SE3Curve_t>(fileName);
  SE3Curve_t se3_from_binary_bezier;
  se3_from_binary_bezier.loadFromBinary<SE3Curve_t>(fileName);
  CompareCurves<SE3Curve_t, SE3Curve_t>(cBezier, se3_from_binary_bezier,
                                        errmsg + " For binary serialization",
                                        error);
}

/**
 * @brief BezierLinearProblemTests test the generation of linear / quadratic
 * problems with variable control points bezier curves
 * @param error
 */

using namespace ndcurves::optimization;

var_pair_t setup_control_points(
    const std::size_t degree, const constraint_flag flag,
    const point3_t& initPos = point3_t(), const point3_t& endPos = point3_t(),
    const constraint_linear& constraints = constraint_linear(3),
    const double totalTime = 1.) {
  problem_definition_t pDef(constraints);
  pDef.init_pos = initPos;
  pDef.end_pos = endPos;
  pDef.flag = flag;
  pDef.totalTime = totalTime;
  pDef.degree = degree;
  problem_data_t pData = setup_control_points<point3_t, double, true>(pDef);
  return std::make_pair(
      pData.variables_,
      std::make_pair(pData.startVariableIndex, pData.numVariables));
}

enum vartype { variable, constant };

bool isVar(const linear_variable_t& var) {
  return !var.isZero() &&
         var.B() == linear_variable_t::matrix_x_t::Identity(3, 3) &&
         var.c() == linear_variable_t::vector_x_t::Zero(3);
}

bool isConstant(const linear_variable_t& var) {
  return var.isZero() ||
         (var.B() == linear_variable_t::matrix_x_t::Zero(3, 3) &&
          var.c() != linear_variable_t::vector_x_t::Zero(3));
}

/*bool isMixed(const linear_variable_t& var)
{
    return var.A_ != linear_variable_t::matrix_t::Zero() &&
           var.b_ != linear_variable_t::point_t::Zero();
}*/

bool checkValue(const linear_variable_t& var, const vartype vart) {
  if (vart == constant)
    return isConstant(var);
  else
    return isVar(var);
}

void checksequence(const T_linear_variable_t& vars, vartype* expected,
                   const std::string testname, bool& error) {
  int i = 0;
  for (CIT_linear_variable_t cit = vars.begin(); cit != vars.end();
       ++cit, ++i) {
    if (!checkValue(*cit, expected[i])) {
      std::cout << "in test: " << testname
                << ": wrong type for variable at position " << i << std::endl;
      error = true;
    }
  }
}

void checkNumVar(const T_linear_variable_t& vars, const std::size_t expected,
                 const std::string testname, bool& error) {
  if (vars.size() != expected) {
    error = true;
    std::cout << "incorrect number of variables in " << testname << "("
              << expected << "," << vars.size() << ")" << std::endl;
  }
}

void checkPair(const pair_size_t pair, const std::size_t start_index,
               const std::size_t num_vars, const std::string testname,
               bool& error) {
  if (pair.first != start_index) {
    error = true;
    std::cout << "incorrect starting index for variablesin "

              << testname << "(" << start_index << "," << pair.first << ")"
              << std::endl;
  }
  if (pair.second != num_vars) {
    error = true;
    std::cout << "incorrect number of identified variablesin " << testname
              << "(" << num_vars << "," << pair.second << ")" << std::endl;
  }
}

void BezierLinearProblemsetup_control_pointsNoConstraint(bool& error) {
  constraint_flag flag = optimization::NONE;
  var_pair_t res_no_constraints = setup_control_points(5, flag);
  T_linear_variable_t& vars = res_no_constraints.first;
  vartype exptecdvars[] = {variable, variable, variable,
                           variable, variable, variable};
  checkNumVar(vars, 6, "setup_control_pointsNoConstraint", error);
  checksequence(vars, exptecdvars, "setup_control_pointsNoConstraint", error);
}

constraint_linear makeConstraint() {
  point3_t init_pos = point3_t(1., 1., 1.);
  constraint_linear cl(3);
  init_pos *= 2;
  cl.init_vel = init_pos;
  init_pos *= 2;
  cl.init_acc = init_pos;
  init_pos *= 2;
  cl.end_acc = init_pos;
  init_pos *= 2;
  cl.end_vel = init_pos;
  return cl;
}

void BezierLinearProblemsetup_control_pointsVarCombinatorialInit(bool& error) {
  constraint_flag flag = optimization::INIT_POS;
  point3_t init_pos = point3_t(1., 1., 1.);
  var_pair_t res = setup_control_points(5, flag, init_pos);
  T_linear_variable_t& vars = res.first;
  vartype exptecdvars[] = {constant, variable, variable,
                           variable, variable, variable};
  checkNumVar(vars, 6, "VarCombinatorialInit", error);
  checksequence(vars, exptecdvars, "VarCombinatorialInit", error);
  checkPair(res.second, 1, 5, "VarCombinatorialInit", error);

  constraint_linear constraints = makeConstraint();
  flag = INIT_POS | INIT_VEL;
  res = setup_control_points(5, flag, init_pos, point3_t(), constraints);
  vars = res.first;
  vartype exptecdvar1[] = {constant, constant, variable,
                           variable, variable, variable};
  checkNumVar(vars, 6, "VarCombinatorialInit", error);
  checksequence(vars, exptecdvar1, "VarCombinatorialInit", error);
  checkPair(res.second, 2, 4, "VarCombinatorialInit", error);

  flag = INIT_POS | INIT_VEL | INIT_ACC;
  res = setup_control_points(5, flag, init_pos, point3_t(), constraints);
  vars = res.first;
  vartype exptecdvar2[] = {constant, constant, constant,
                           variable, variable, variable};
  checkNumVar(vars, 6, "VarCombinatorialInit", error);
  checksequence(vars, exptecdvar2, "VarCombinatorialInit", error);
  checkPair(res.second, 3, 3, "VarCombinatorialInit", error);

  flag = INIT_VEL;
  res = setup_control_points(5, flag, init_pos, point3_t(), constraints);
  vars = res.first;
  vartype exptecdvar3[] = {variable, variable, variable,
                           variable, variable, variable};
  checkNumVar(vars, 6, "VarCombinatorialInit", error);
  checksequence(vars, exptecdvar3, "VarCombinatorialInit", error);
  checkPair(res.second, 0, 6, "VarCombinatorialInit", error);

  flag = INIT_ACC;
  res = setup_control_points(5, flag, init_pos, point3_t(), constraints);
  vars = res.first;
  vartype exptecdvar4[] = {variable, variable, variable,
                           variable, variable, variable};
  checkNumVar(vars, 6, "VarCombinatorialInit", error);
  checksequence(vars, exptecdvar4, "VarCombinatorialInit", error);
  checkPair(res.second, 0, 6, "VarCombinatorialInit", error);

  flag = INIT_ACC | INIT_VEL;
  res = setup_control_points(5, flag, init_pos, point3_t(), constraints);
  vars = res.first;
  vartype exptecdvar5[] = {variable, variable, variable,
                           variable, variable, variable};
  checkNumVar(vars, 6, "VarCombinatorialInit", error);
  checksequence(vars, exptecdvar5, "VarCombinatorialInit", error);
  checkPair(res.second, 0, 6, "VarCombinatorialInit", error);

  bool err = true;
  try {
    flag = INIT_POS | INIT_VEL;
    res = setup_control_points(1, flag, init_pos, point3_t(), constraints);
  } catch (...) {
    err = false;
  }
  if (err) {
    error = true;
    std::cout << "exception should be raised when degree of bezier curve is "
                 "not high enough to handle constraints "
              << std::endl;
  }
}

void BezierLinearProblemsetup_control_pointsVarCombinatorialEnd(bool& error) {
  constraint_flag flag = optimization::END_POS;
  point3_t init_pos = point3_t(1., 1., 1.);
  var_pair_t res = setup_control_points(5, flag, init_pos, init_pos);
  T_linear_variable_t& vars = res.first;
  vartype exptecdvars[] = {variable, variable, variable,
                           variable, variable, constant};
  checkNumVar(vars, 6, "VarCombinatorialEnd", error);
  checksequence(vars, exptecdvars, "VarCombinatorialEnd", error);
  checkPair(res.second, 0, 5, "VarCombinatorialEnd", error);

  constraint_linear constraints = makeConstraint();
  flag = END_POS | END_VEL;
  res = setup_control_points(5, flag, init_pos, init_pos, constraints);
  vars = res.first;
  vartype exptecdvar1[] = {variable, variable, variable,
                           variable, constant, constant};
  checkNumVar(vars, 6, "VarCombinatorialEnd", error);
  checksequence(vars, exptecdvar1, "VarCombinatorialEnd", error);
  checkPair(res.second, 0, 4, "VarCombinatorialEnd", error);

  flag = END_POS | END_VEL | END_ACC;
  res = setup_control_points(5, flag, init_pos, init_pos, constraints);
  vars = res.first;
  vartype exptecdvar2[] = {variable, variable, variable,
                           constant, constant, constant};
  checkNumVar(vars, 6, "VarCombinatorialEnd", error);
  checksequence(vars, exptecdvar2, "VarCombinatorialEnd", error);
  checkPair(res.second, 0, 3, "VarCombinatorialEnd", error);

  flag = END_VEL;
  res = setup_control_points(5, flag, init_pos, init_pos, constraints);
  vars = res.first;
  vartype exptecdvar3[] = {variable, variable, variable,
                           variable, variable, variable};
  checkNumVar(vars, 6, "VarCombinatorialEnd", error);
  checksequence(vars, exptecdvar3, "VarCombinatorialEnd", error);
  checkPair(res.second, 0, 6, "VarCombinatorialEnd", error);

  flag = END_ACC;
  res = setup_control_points(5, flag, init_pos, init_pos, constraints);
  vars = res.first;
  vartype exptecdvar4[] = {variable, variable, variable,
                           variable, variable, variable};
  checkNumVar(vars, 6, "VarCombinatorialEnd", error);
  checksequence(vars, exptecdvar4, "VarCombinatorialEnd", error);
  checkPair(res.second, 0, 6, "VarCombinatorialEnd", error);

  flag = END_ACC | END_VEL;
  res = setup_control_points(5, flag, init_pos, init_pos, constraints);
  vars = res.first;
  vartype exptecdvar5[] = {variable, variable, variable,
                           variable, variable, variable};
  checkNumVar(vars, 6, "VarCombinatorialEnd", error);
  checksequence(vars, exptecdvar5, "VarCombinatorialEnd", error);
  checkPair(res.second, 0, 6, "VarCombinatorialEnd", error);

  bool err = true;
  try {
    flag = END_ACC | END_VEL;
    res = setup_control_points(1, flag, init_pos, point3_t(), constraints);
  } catch (...) {
    err = false;
  }
  if (err) {
    error = true;
    std::cout << "exception should be raised when degree of bezier curve is "
                 "not high enough to handle constraints "
              << std::endl;
  }
}

void BezierLinearProblemsetup_control_pointsVarCombinatorialMix(bool& error) {
  constraint_flag flag = END_POS | INIT_POS;
  point3_t init_pos = point3_t(1., 1., 1.);
  var_pair_t res = setup_control_points(5, flag, init_pos, init_pos);
  T_linear_variable_t& vars = res.first;
  vartype exptecdvars[] = {constant, variable, variable,
                           variable, variable, constant};
  checkNumVar(vars, 6, "VarCombinatorialMix", error);
  checksequence(vars, exptecdvars, "VarCombinatorialMix", error);
  checkPair(res.second, 1, 4, "VarCombinatorialMix", error);

  constraint_linear constraints = makeConstraint();
  flag = END_POS | END_VEL | INIT_VEL | INIT_POS;
  res = setup_control_points(5, flag, init_pos, init_pos, constraints);
  vars = res.first;
  vartype exptecdvar1[] = {constant, constant, variable,
                           variable, constant, constant};
  checkNumVar(vars, 6, "VarCombinatorialMix", error);
  checksequence(vars, exptecdvar1, "VarCombinatorialMix", error);
  checkPair(res.second, 2, 2, "VarCombinatorialMix", error);

  flag = END_POS | END_VEL | END_ACC | INIT_VEL | INIT_POS;
  res = setup_control_points(5, flag, init_pos, init_pos, constraints);
  vars = res.first;
  vartype exptecdvar2[] = {constant, constant, variable,
                           constant, constant, constant};
  checkNumVar(vars, 6, "VarCombinatorialMix", error);
  checksequence(vars, exptecdvar2, "VarCombinatorialMix", error);
  checkPair(res.second, 2, 1, "VarCombinatorialMix", error);

  flag = ALL;
  res = setup_control_points(8, flag, init_pos, init_pos, constraints);
  vars = res.first;
  vartype exptecdvar3[] = {constant, constant, constant, constant, variable,
                           constant, constant, constant, constant};
  checkNumVar(vars, 9, "VarCombinatorialMix", error);
  checksequence(vars, exptecdvar3, "VarCombinatorialMix", error);
  checkPair(res.second, 4, 1, "VarCombinatorialMix", error);

  flag = END_VEL | END_ACC | INIT_VEL;
  res = setup_control_points(5, flag, init_pos, init_pos, constraints);
  vars = res.first;
  vartype exptecdvar4[] = {variable, variable, variable,
                           variable, variable, variable};
  checkNumVar(vars, 6, "VarCombinatorialMix", error);
  checksequence(vars, exptecdvar4, "VarCombinatorialMix", error);
  checkPair(res.second, 0, 6, "VarCombinatorialMix", error);

  flag = END_VEL | INIT_VEL;
  res = setup_control_points(5, flag, init_pos, init_pos, constraints);
  vars = res.first;
  vartype exptecdvar5[] = {variable, variable, variable,
                           variable, variable, variable};
  checkNumVar(vars, 6, "VarCombinatorialMix", error);
  checksequence(vars, exptecdvar5, "VarCombinatorialMix", error);
  checkPair(res.second, 0, 6, "VarCombinatorialMix", error);

  bool err = true;
  try {
    flag = ALL;
    res = setup_control_points(5, flag, init_pos, init_pos, constraints);
  } catch (...) {
    err = false;
  }
  if (err) {
    error = true;
    std::cout << "exception should be raised when degree of bezier curve is "
                 "not high enough to handle constraints "
              << std::endl;
  }
}

void BezierLinearProblemInitInequalities(bool& error) {
  constraint_flag flag = INIT_POS | END_POS;
  point3_t init_pos = point3_t(1., 1., 1.);
  var_pair_t res = setup_control_points(5, flag, init_pos, init_pos);
  T_linear_variable_t& vars = res.first;
  vartype exptecdvars[] = {constant, variable, variable,
                           variable, variable, constant};
  checkNumVar(vars, 6, "VarCombinatorialMix", error);
  checksequence(vars, exptecdvars, "VarCombinatorialMix", error);
  checkPair(res.second, 1, 4, "VarCombinatorialMix", error);

  constraint_linear constraints = makeConstraint();
  flag = END_POS | END_VEL | INIT_VEL | INIT_POS;
  res = setup_control_points(5, flag, init_pos, init_pos, constraints);
  vars = res.first;
  vartype exptecdvar1[] = {constant, constant, variable,
                           variable, constant, constant};
  checkNumVar(vars, 6, "VarCombinatorialMix", error);
  checksequence(vars, exptecdvar1, "VarCombinatorialMix", error);
  checkPair(res.second, 2, 2, "VarCombinatorialMix", error);

  flag = END_POS | END_VEL | END_ACC | INIT_VEL | INIT_POS;
  res = setup_control_points(5, flag, init_pos, init_pos, constraints);
  vars = res.first;
  vartype exptecdvar2[] = {constant, constant, variable,
                           constant, constant, constant};
  checkNumVar(vars, 6, "VarCombinatorialMix", error);
  checksequence(vars, exptecdvar2, "VarCombinatorialMix", error);
  checkPair(res.second, 2, 1, "VarCombinatorialMix", error);

  flag = ALL;
  res = setup_control_points(6, flag, init_pos, init_pos, constraints);
  vars = res.first;
  vartype exptecdvar3[] = {constant, constant, constant, variable,
                           constant, constant, constant};
  checkNumVar(vars, 7, "VarCombinatorialMix", error);
  checksequence(vars, exptecdvar3, "VarCombinatorialMix", error);
  checkPair(res.second, 3, 1, "VarCombinatorialMix", error);

  flag = END_VEL | END_ACC | INIT_VEL;
  res = setup_control_points(5, flag, init_pos, init_pos, constraints);
  vars = res.first;
  vartype exptecdvar4[] = {variable, variable, variable,
                           variable, variable, variable};
  checkNumVar(vars, 6, "VarCombinatorialMix", error);
  checksequence(vars, exptecdvar4, "VarCombinatorialMix", error);
  checkPair(res.second, 0, 6, "VarCombinatorialMix", error);

  flag = END_VEL | INIT_VEL;
  res = setup_control_points(5, flag, init_pos, init_pos, constraints);
  vars = res.first;
  vartype exptecdvar5[] = {variable, variable, variable,
                           variable, variable, variable};
  checkNumVar(vars, 6, "VarCombinatorialMix", error);
  checksequence(vars, exptecdvar5, "VarCombinatorialMix", error);
  checkPair(res.second, 0, 6, "VarCombinatorialMix", error);

  bool err = true;
  try {
    flag = ALL;
    res = setup_control_points(5, flag, init_pos, init_pos, constraints);
  } catch (...) {
    err = false;
  }
  if (err) {
    error = true;
    std::cout << "exception should be raised when degree of bezier curve is "
                 "not high enough to handle constraints "
              << std::endl;
  }
}

void BezierLinearProblemsetupLoadProblem(bool& /*error*/) {
  problem_definition_t pDef = loadproblem(TEST_DATA_PATH "test.pb");
  // problem_data_t pData = setup_control_points<point_t, 3, double>(pDef);
  generate_problem<point3_t, double, true>(pDef, VELOCITY);
  // initInequalityMatrix<point_t,3,double>(pDef,pData,prob);
}

void testOperatorEqual(bool& error) {
  // test with a C2 polynomial :
  pointX_t zeros = point3_t(0., 0., 0.);
  pointX_t p0 = point3_t(0., 1., 0.);
  pointX_t p1 = point3_t(1., 2., -3.);
  pointX_t dp0 = point3_t(-8., 4., 6.);
  pointX_t dp1 = point3_t(10., -10., 10.);
  pointX_t ddp0 = point3_t(-1., 7., 4.);
  pointX_t ddp1 = point3_t(12., -8., 2.5);
  double min = 0.5;
  double max = 2.;
  polynomial_t polC2_1 = polynomial_t(p0, dp0, ddp0, p1, dp1, ddp1, min, max);
  polynomial_t polC2_2 = polynomial_t(p0, dp0, ddp0, p1, dp1, ddp1, min, max);
  polynomial_t polC2_3(polC2_1);
  // std::cout<<"Should call polynomial method : "<<std::endl;
  if (polC2_1 != polC2_2) {
    std::cout << "polC2_1 and polC2_2 should be equals" << std::endl;
    error = true;
  }
  if (polC2_1 != polC2_3) {
    std::cout << "polC2_1 and polC2_3 should be equals" << std::endl;
    error = true;
  }

  polynomial_t polC2_4 = polynomial_t(p1, dp0, ddp0, p1, dp1, ddp0, min, max);
  if (polC2_1 == polC2_4) {
    std::cout << "polC2_1 and polC2_4 should not be equals" << std::endl;
    error = true;
  }

  // test with bezier
  point3_t a(1, 2, 3);
  point3_t b(2, 3, 4);
  point3_t c(3, 4, 5);
  point3_t d(3, 6, 7);
  point3_t e(3, 61, 7);
  point3_t f(3, 56, 7);
  point3_t g(3, 36, 7);
  point3_t h(43, 6, 7);
  point3_t i(3, 6, 77);
  std::vector<point3_t> control_points;
  control_points.push_back(a);
  control_points.push_back(b);
  control_points.push_back(c);
  control_points.push_back(d);
  control_points.push_back(e);
  control_points.push_back(f);
  control_points.push_back(g);
  control_points.push_back(h);
  control_points.push_back(i);
  bezier_t::num_t T_min = 1.0;
  bezier_t::num_t T_max = 3.0;
  bezier_t bc_0(control_points.begin(), control_points.end(), T_min, T_max);
  bezier3_t bc3_0(control_points.begin(), control_points.end(), T_min, T_max);
  bezier_t bc_1(bc_0);
  // std::cout<<"Should call Bezier method : "<<std::endl;
  if (bc_1 != bc_0) {
    std::cout << "bc_0 and bc_1 should be equals" << std::endl;
    error = true;
  }
  std::vector<point3_t> control_points2;
  control_points2.push_back(a);
  control_points2.push_back(b);
  control_points2.push_back(c);
  control_points2.push_back(d);
  bezier_t bc_2(control_points2.begin(), control_points2.end(), T_min, T_max);
  bezier3_t bc3_2(control_points2.begin(), control_points2.end(), T_min, T_max);
  if (bc_2 == bc_0) {
    std::cout << "bc_0 and bc_2 should not be equals" << std::endl;
    error = true;
  }

  point3_t e3(3, 61.9, 7);
  point3_t g3(-3, 36, 7);
  std::vector<point3_t> control_points3;
  control_points3.push_back(a);
  control_points3.push_back(b);
  control_points3.push_back(c);
  control_points3.push_back(d);
  control_points3.push_back(e3);
  control_points3.push_back(f);
  control_points3.push_back(g3);
  control_points3.push_back(h);
  control_points3.push_back(i);
  bezier_t bc_0_3(control_points3.begin(), control_points3.end(), T_min, T_max);
  if (bc_0_3 == bc_0) {
    std::cout << "bc_0_3 and bc_0 should not be equals" << std::endl;
    error = true;
  }
  polynomial_t pol_2 = polynomial_from_curve<polynomial_t>(bc_2);
  polynomial3_t pol3_2 = polynomial_from_curve<polynomial3_t>(bc3_2);
  bezier_t bc_3 = bezier_from_curve<bezier_t>(pol_2);
  if (bc_2 != bc_3) {
    std::cout << "bc_2 and bc_3 should be equals" << std::endl;
    error = true;
  }

  // test bezier / polynomial
  polynomial_t pol_0 = polynomial_from_curve<polynomial_t>(bc_0);
  CompareCurves<polynomial_t, bezier_t>(pol_0, bc_0, "compare pol_0 and bc_0",
                                        error);
  // std::cout<<"Should call curve_abc method : "<<std::endl;
  if (!bc_0.isEquivalent(&pol_0)) {
    std::cout << "bc_0 and pol_0 should be equivalent" << std::endl;
    error = true;
  }

  // test with hermite :
  point3_t ch_p0(1, 2, 3);
  point3_t ch_m0(2, 3, 4);
  point3_t ch_p1(3, 4, 5);
  point3_t ch_m1(3, 6, 7);
  pair_point_tangent_t pair0(ch_p0, ch_m0);
  pair_point_tangent_t pair1(ch_p1, ch_m1);
  t_pair_point_tangent_t ch_control_points;
  ch_control_points.push_back(pair0);
  ch_control_points.push_back(pair1);
  std::vector<double> time_control_points;
  time_control_points.push_back(T_min);
  time_control_points.push_back(T_max);
  cubic_hermite_spline_t chs0(ch_control_points.begin(),
                              ch_control_points.end(), time_control_points);
  cubic_hermite_spline_t chs1(ch_control_points.begin(),
                              ch_control_points.end(), time_control_points);
  cubic_hermite_spline_t chs2(chs0);
  point3_t ch_p2(3.1, 4, 5);
  point3_t ch_m2(3, 6.5, 6.);
  pair_point_tangent_t pair2(ch_p2, ch_m2);
  t_pair_point_tangent_t ch_control_points2;
  ch_control_points2.push_back(pair0);
  ch_control_points2.push_back(pair2);
  cubic_hermite_spline_t chs3(ch_control_points2.begin(),
                              ch_control_points2.end(), time_control_points);
  // std::cout<<"Should call hermite method : "<<std::endl;
  if (chs0 != chs1) {
    std::cout << "chs0 and chs1 should be equals" << std::endl;
    error = true;
  }
  if (chs0 != chs2) {
    std::cout << "chs0 and chs2 should be equals" << std::endl;
    error = true;
  }
  if (chs0 == chs3) {
    std::cout << "chs0 and chs3 should not be equals" << std::endl;
    error = true;
  }

  //  // test bezier / hermite
  bezier_t bc_ch = bezier_from_curve<bezier_t>(chs0);
  // std::cout<<"Should call curve_abc method : "<<std::endl;
  if (!chs0.isEquivalent(&bc_ch)) {
    std::cout << "chs0 and bc_ch should be equivalent" << std::endl;
    error = true;
  }
  // test polynomial / hermite
  polynomial_t pol_ch = polynomial_from_curve<polynomial_t>(chs0);
  if (!chs0.isEquivalent(&pol_ch)) {
    std::cout << "chs0 and pol_ch should be equivalent" << std::endl;
    error = true;
  }

  // SO3
  quaternion_t q0(1, 0, 0, 0);
  quaternion_t q1(0.7071, 0.7071, 0, 0);
  q0.normalize();
  q1.normalize();
  const double tMin = 0.;
  const double tMax = 1.5;
  SO3Linear_t so3Traj1(q0, q1, tMin, tMax);
  SO3Linear_t so3Traj2(q0, q1, tMin, tMax);
  SO3Linear_t so3TrajMatrix1(q0.toRotationMatrix(), q1.toRotationMatrix(), tMin,
                             tMax);
  SO3Linear_t so3TrajMatrix2(so3TrajMatrix1);
  quaternion_t q2(0.7071, 0., 0.7071, 0);
  q2.normalize();
  SO3Linear_t so3Traj3(q0, q2, tMin, tMax);
  // std::cout<<"Should call SO3 method : "<<std::endl;
  if (so3Traj1 != so3Traj2) {
    std::cout << "so3Traj1 and so3Traj2 should be equals" << std::endl;
    error = true;
  }
  if (so3Traj1 != so3TrajMatrix1) {
    std::cout << "so3Traj1 and so3TrajMatrix1 should be equals" << std::endl;
    error = true;
  }
  if (so3Traj1 != so3TrajMatrix2) {
    std::cout << "so3Traj1 and so3TrajMatrix2 should be equals" << std::endl;
    error = true;
  }
  if (so3Traj1 == so3Traj3) {
    std::cout << "so3Traj1 and so3Traj3 should not be equals" << std::endl;
    error = true;
  }

  // test from pointer :
  curve_ptr_t c_ptr1(new bezier_t(bc_0));
  curve_abc_t* c_ptr2 = new bezier_t(bc_0);
  curve_ptr_t c_ptr3(new polynomial_t(pol_0));
  // std::cout<<"Should call bezier method : "<<std::endl;
  if (!c_ptr1->isApprox(c_ptr2)) {
    std::cout << "c_ptr1 and c_ptr2 should be approx" << std::endl;
    error = true;
  }
  // std::cout<<"Should call curve_abc method : "<<std::endl;
  if (!c_ptr2->isEquivalent(c_ptr3.get())) {
    std::cout << "c_ptr2 and c_ptr3 should be equivalent" << std::endl;
    error = true;
  }
  if (c_ptr1->isApprox(c_ptr3.get())) {
    std::cout << "c_ptr1 and c_ptr3 should not be approx" << std::endl;
    error = true;
  }

  // SE3
  polynomial3_t pol3_0 = polynomial_from_curve<polynomial3_t>(bc3_0);
  std::shared_ptr<bezier3_t> translation_bezier =
      std::make_shared<bezier3_t>(bc3_0);
  std::shared_ptr<bezier3_t> translation_bezier2 =
      std::make_shared<bezier3_t>(bc3_0);
  std::shared_ptr<polynomial3_t> translation_polynomial =
      std::make_shared<polynomial3_t>(pol3_0);
  SE3Curve_t se3_bezier1 = SE3Curve_t(translation_bezier, q0.toRotationMatrix(),
                                      q1.toRotationMatrix());
  SE3Curve_t se3_bezier12 = SE3Curve_t(
      translation_bezier2, q0.toRotationMatrix(), q1.toRotationMatrix());
  SE3Curve_t se3_pol1 = SE3Curve_t(
      translation_polynomial, q0.toRotationMatrix(), q1.toRotationMatrix());
  SE3Curve_t se3_bezier2(se3_bezier1);
  SE3Curve_t se3_bezier3 = SE3Curve_t(translation_bezier, q0.toRotationMatrix(),
                                      q2.toRotationMatrix());
  std::shared_ptr<polynomial3_t> translation_polynomial2 =
      std::make_shared<polynomial3_t>(pol3_2);
  SE3Curve_t se3_pol2 = SE3Curve_t(
      translation_polynomial2, q0.toRotationMatrix(), q1.toRotationMatrix());
  // std::cout<<"Should call se3 method : "<<std::endl;
  if (se3_bezier1 == se3_pol1) {
    std::cout << "se3_bezier1 and se3_pol1 should not be equals" << std::endl;
    error = true;
  }
  if (se3_bezier1.isApprox(se3_pol1)) {
    std::cout << "se3_bezier1 and se3_pol1 should not be approx" << std::endl;
    error = true;
  }
  // std::cout<<"Should call curve_abc : "<<std::endl;
  if (!se3_bezier1.isEquivalent(&se3_pol1)) {
    std::cout << "se3_bezier1 and se3_pol1 should be equivalent" << std::endl;
    error = true;
  }
  // std::cout<<"Should call se3 method : "<<std::endl;
  if (se3_bezier1 != se3_bezier2) {
    std::cout << "se3_bezier1 and se3_bezier2 should be equals" << std::endl;
    error = true;
  }
  // std::cout<<"Should call se3 -> bezier / SO3 method : "<<std::endl;
  if (se3_bezier1 != se3_bezier12) {
    std::cout << "se3_bezier1 and se3_bezier12 should be equals" << std::endl;
    error = true;
  }
  // std::cout<<"Should call se3 -> curve_abc : "<<std::endl;
  if (se3_bezier1 == se3_pol2) {
    std::cout << "se3_bezier1 and se3_pol2 should not be equals" << std::endl;
    error = true;
  }
  // std::cout<<"Should call se3 -> so3 method : "<<std::endl;
  if (se3_bezier1 == se3_bezier3) {
    std::cout << "se3_bezier1 and se3_bezier3 should not be equals"
              << std::endl;
    error = true;
  }

  // Piecewises
  point3_t a0(1, 2, 3);
  point3_t b0(2, 3, 4);
  point3_t c0(3, 4, 5);
  point3_t d0(4, 5, 6);
  std::vector<point3_t> params0;
  std::vector<point3_t> params1;
  params0.push_back(a0);  // bezier between [0,1]
  params0.push_back(b0);
  params0.push_back(c0);
  params0.push_back(d0);
  params1.push_back(d0);  // bezier between [1,2]
  params1.push_back(c0);
  params1.push_back(b0);
  params1.push_back(a0);
  std::shared_ptr<bezier_t> bc0_ptr =
      std::make_shared<bezier_t>(params0.begin(), params0.end(), 0., 1.);
  std::shared_ptr<bezier_t> bc1_ptr =
      std::make_shared<bezier_t>(params1.begin(), params1.end(), 1., 2.);
  piecewise_t pc_C0(bc0_ptr);
  pc_C0.add_curve_ptr(bc1_ptr);
  piecewise_t pc_C1(bc0_ptr);
  pc_C1.add_curve_ptr(bc1_ptr);
  piecewise_t pc_C2(pc_C0);
  piecewise_t pc_C3(bc0_ptr);
  piecewise_t pc_C4(bc0_ptr);
  std::shared_ptr<bezier_t> bc2_ptr =
      std::make_shared<bezier_t>(params0.begin(), params0.end(), 1., 2.);
  pc_C4.add_curve_ptr(bc2_ptr);
  // std::cout<<"Should call piecewise method -> bezier , bezier: "<<std::endl;
  if (pc_C0 != pc_C1) {
    std::cout << "pc_C0 and pc_C1 should be equals" << std::endl;
    error = true;
  }
  // std::cout<<"Should call piecewise method -> bezier , bezier: "<<std::endl;
  if (pc_C0 != pc_C2) {
    std::cout << "pc_C0 and pc_C2 should be equals" << std::endl;
    error = true;
  }
  // std::cout<<"Should call piecewise method: "<<std::endl;
  if (pc_C0 == pc_C3) {
    std::cout << "pc_C0 and pc_C3 should not be equals" << std::endl;
    error = true;
  }
  // std::cout<<"Should call piecewise method -> bezier , bezier: "<<std::endl;
  if (pc_C0 == pc_C4) {
    std::cout << "pc_C0 and pc_C4 should not be equals" << std::endl;
    error = true;
  }
  // piecewise with mixed curves types
  // std::cout<<"Should call piecewise method: "<<std::endl;
  piecewise_t pc_C5 =
      pc_C0.convert_piecewise_curve_to_polynomial<polynomial_t>();
  if (pc_C0 == pc_C5) {
    std::cout << "pc_C0 and pc_C5 should be not equals" << std::endl;
    error = true;
  }
  // std::cout<<"Should call curve_abc method: "<<std::endl;
  if (!pc_C0.isEquivalent(&pc_C5)) {
    std::cout << "pc_C0 and pc_C5 should be equivalents" << std::endl;
    error = true;
  }

  // piecewise se3 :
  piecewise_SE3_t pc_se3_1;
  pc_se3_1.add_curve(se3_pol1);
  point3_t p_init_se3(
      translation_polynomial->operator()(translation_polynomial->max()));
  point3_t dp_init_se3(
      translation_polynomial->derivate(translation_polynomial->max(), 1));
  point3_t p_end_se3(1, -2, 6);
  point3_t dp_end_se3(3.5, 2.5, -9);
  std::shared_ptr<polynomial3_t> translation_pol3 =
      std::make_shared<polynomial3_t>(p_init_se3, dp_init_se3, p_end_se3,
                                      dp_end_se3, translation_polynomial->max(),
                                      translation_polynomial->max() + 2.5);
  curve_SE3_ptr_t se3_pol_3(new SE3Curve_t(
      translation_pol3, q1.toRotationMatrix(), q2.toRotationMatrix()));
  pc_se3_1.add_curve_ptr(se3_pol_3);
  piecewise_SE3_t pc_se3_2(pc_se3_1);
  piecewise_SE3_t pc_se3_3(std::make_shared<SE3Curve_t>(se3_pol1));
  pc_se3_3.add_curve_ptr(se3_pol_3);
  piecewise_SE3_t pc_se3_4(std::make_shared<SE3Curve_t>(se3_pol2));
  pc_se3_4.add_curve_ptr(se3_pol_3);
  // std::cout<<"Should call piecewise method -> SE3 , SE3: "<<std::endl;
  if (pc_se3_1 != pc_se3_2) {
    std::cout << "pc_se3_1 and pc_se3_2 should be equals" << std::endl;
    error = true;
  }
  // std::cout<<"Should call piecewise method -> SE3 , SE3: "<<std::endl;
  if (pc_se3_1 != pc_se3_3) {
    std::cout << "pc_se3_1 and pc_se3_3 should be equals" << std::endl;
    error = true;
  }
  // std::cout<<"Should call piecewise method -> SE3  -> polynomial :
  // "<<std::endl;
  if (pc_se3_1 == pc_se3_4) {
    std::cout << "pc_se3_1 and pc_se3_3 should not be equals" << std::endl;
    error = true;
  }
}

int main(int /*argc*/, char** /*argv[]*/) {
  std::cout << "performing tests... \n";
  bool error = false;
  PolynomialCubicFunctionTest(error);
  ExactCubicNoErrorTest(error);
  ExactCubicPointsCrossedTest(
      error);  // checks that given wayPoints are crossed
  ExactCubicTwoPointsTest(error);
  ExactCubicOneDimTest(error);
  ExactCubicVelocityConstraintsTest(error);
  EffectorTrajectoryTest(error);
  EffectorSplineRotationNoRotationTest(error);
  EffectorSplineRotationRotationTest(error);
  TestReparametrization(error);
  EffectorSplineRotationWayPointRotationTest(error);
  BezierCurveTest(error);
  BezierDerivativeCurveTest(error);
  BezierDerivativeCurveConstraintTest(error);
  BezierCurveTestCompareHornerAndBernstein(error);
  BezierDerivativeCurveTimeReparametrizationTest(error);
  BezierEvalDeCasteljau(error);
  BezierSplitCurve(error);
  BezierElevate(error);
  CubicHermitePairsPositionDerivativeTest(error);
  piecewiseCurveTest(error);
  PiecewisePolynomialCurveFromDiscretePoints(error);
  PiecewisePolynomialCurveFromFile(error);
  toPolynomialConversionTest(error);
  cubicConversionTest(error);
  curveAbcDimDynamicTest(error);
  serializationCurvesTest(error);
  polynomialFromBoundaryConditions(error);
  so3LinearTest(error);
  SO3serializationTest(error);
  se3CurveTest(error);
  Se3serializationTest(error);
  BezierLinearProblemsetup_control_pointsNoConstraint(error);
  BezierLinearProblemsetup_control_pointsVarCombinatorialInit(error);
  BezierLinearProblemsetup_control_pointsVarCombinatorialEnd(error);
  BezierLinearProblemsetup_control_pointsVarCombinatorialMix(error);
  BezierLinearProblemsetupLoadProblem(error);
  testOperatorEqual(error);

  if (error) {
    std::cout << "There were some errors\n";
    return -1;
  } else {
    std::cout << "no errors found \n";
    return 0;
  }
}
