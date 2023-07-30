#define BOOST_TEST_MODULE test_operations

#include <boost/test/included/unit_test.hpp>

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

using namespace ndcurves;

BOOST_AUTO_TEST_SUITE(BOOST_TEST_MODULE)

void compDouble(const double a, const double b) {
#if BOOST_VERSION <= 105800
  BOOST_CHECK(std::abs(a - b) <= 0.001);
#else
  BOOST_TEST(a == b, boost::test_tools::tolerance(0.001));
#endif
}

BOOST_AUTO_TEST_CASE(crossPoductBezier) {
  t_pointX_t vec1;
  t_pointX_t vec2;
  for (int i = 0; i < 4; ++i) {
    vec1.push_back(Eigen::Vector3d::Random());
    vec2.push_back(Eigen::Vector3d::Random());
  }
  for (int i = 0; i < 3; ++i) {
    vec1.push_back(Eigen::Vector3d::Random());
  }
  bezier_t p1(vec1.begin(), vec1.end(), 0., 1.);
  bezier_t p2(vec2.begin(), vec2.end(), 0., 1.);

  bezier_t pCross(p1.cross(p2));
  for (double i = 0.; i <= 100.; ++i) {
    double dt = i / 100.;
    Eigen::Vector3d v1 = p1(dt);
    Eigen::Vector3d v2 = p2(dt);
    compDouble((pCross(dt) - v1.cross(v2)).norm(), 0.);
  }
}

BOOST_AUTO_TEST_CASE(bezierOperations) {
  t_pointX_t vec1;
  t_pointX_t vec2;
  for (int i = 0; i < 3; ++i) {
    vec1.push_back(Eigen::Vector3d::Random());
    vec2.push_back(Eigen::Vector3d::Random());
  }
  for (int i = 0; i < 4; ++i) {
    vec1.push_back(Eigen::Vector3d::Random());
  }
  bezier_t p1(vec1.begin(), vec1.end(), 0.2, .8);
  bezier_t p2(vec2.begin(), vec2.end(), 0.2, .8);
  bezier_t p1Dev = p1.compute_derivate(2);
  bezier_t p3(vec2.begin(), vec2.end(), 0.2, 0.5);
  bezier_t p4(vec2.begin(), vec2.end(), 0.1, .8);
  bezier_t p5(vec2.begin(), vec2.end(), 0.1, .5);
  double k = 10.2;

  BOOST_CHECK_THROW(p1 + p3, std::exception);
  BOOST_CHECK_THROW(p1 - p3, std::exception);
  BOOST_CHECK_THROW(p1 + p4, std::exception);
  BOOST_CHECK_THROW(p1 - p4, std::exception);
  BOOST_CHECK_THROW(p1 + p5, std::exception);
  BOOST_CHECK_THROW(p1 - p5, std::exception);

  bezier_t pSum = p1 + p2;
  bezier_t pSumR = p2 + p1;
  bezier_t pSumDev = p1Dev + p2;
  bezier_t pSumRDev = p2 + p1Dev;
  bezier_t pSub = p1 - p2;
  bezier_t pSubR = p2 - p1;
  bezier_t pSubDev = p1Dev - p2;
  bezier_t pSubRDev = p2 - p1Dev;
  bezier_t pdiv = p1 / k;
  bezier_t pMul = p1 * k;
  bezier_t pMulR = k * p1;
  bezier_t pNeg = -p1;
  for (double i = 20; i <= 80.; ++i) {
    double dt = i / 100.;
    compDouble((pSum(dt) - (p1(dt) + p2(dt))).norm(), 0.);
    compDouble((pSumR(dt) - (p1(dt) + p2(dt))).norm(), 0.);
    compDouble((pSumDev(dt) - (p1Dev(dt) + p2(dt))).norm(), 0.);
    compDouble((pSumRDev(dt) - (p1Dev(dt) + p2(dt))).norm(), 0.);
    compDouble((pSub(dt) - (p1(dt) - p2(dt))).norm(), 0.);
    compDouble((pSubR(dt) - (p2(dt) - p1(dt))).norm(), 0.);
    compDouble((pSubDev(dt) - (p1Dev(dt) - p2(dt))).norm(), 0.);
    compDouble((pSubRDev(dt) - (p2(dt) - p1Dev(dt))).norm(), 0.);
    compDouble((pMul(dt) - p1(dt) * k).norm(), 0.);
    compDouble((pMulR(dt) - p1(dt) * k).norm(), 0.);
    compDouble((pdiv(dt) - p1(dt) / k).norm(), 0.);
    compDouble((pNeg(dt) + p1(dt)).norm(), 0);
  }

  pSum = bezier_t(p1);
  pSum += p2;
  pSub = p1;
  pSub -= p2;
  pdiv = p1;
  pdiv /= k;
  pMul = p1;
  pMul *= k;
  for (double i = 20; i <= 80.; ++i) {
    double dt = i / 100.;
    compDouble((pSum(dt) - (p1(dt) + p2(dt))).norm(), 0.);
    compDouble((pSub(dt) - (p1(dt) - p2(dt))).norm(), 0.);
    compDouble((pMul(dt) - p1(dt) * k).norm(), 0.);
    compDouble((pdiv(dt) - p1(dt) / k).norm(), 0.);
  }
}

BOOST_AUTO_TEST_CASE(bezierPointOperations) {
  t_pointX_t vec1;
  for (int i = 0; i < 6; ++i) {
    vec1.push_back(Eigen::Vector3d::Random());
  }

  bezier_t p1(vec1.begin(), vec1.end(), 0., 1.);
  Eigen::Vector3d point = bezier_t::point_t::Random(3);

  bezier_t pSum = p1 + point;
  bezier_t pSumR = point + p1;
  bezier_t pSub = p1 - point;
  bezier_t pSubR = point - p1;
  bezier_t pcross = p1.cross(point);
  for (double i = 0.; i <= 100.; ++i) {
    double dt = i / 100.;
    compDouble((pSum(dt) - (p1(dt) + point)).norm(), 0.);
    compDouble((pSumR(dt) - (p1(dt) + point)).norm(), 0.);
    compDouble((pSub(dt) - (p1(dt) - point)).norm(), 0.);
    compDouble((pSubR(dt) - (point - p1(dt))).norm(), 0.);
    Eigen::Vector3d p1dt = p1(dt);
    compDouble((pcross(dt) - (p1dt.cross(point))).norm(), 0.);
  }
}

BOOST_AUTO_TEST_CASE(crossPoductLinearVariable) {
  linear_variable_t l1(Eigen::Matrix3d::Identity() * 5.,
                       Eigen::Vector3d::Random());
  linear_variable_t l2(Eigen::Matrix3d::Identity() * 1.,
                       Eigen::Vector3d::Random());
  linear_variable_t lE(Eigen::Matrix3d::Random(), Eigen::Vector3d::Random());
  BOOST_CHECK_THROW(l1.cross(lE), std::exception);
  BOOST_CHECK_THROW(lE.cross(l1), std::exception);
  linear_variable_t lcross = l1.cross(l2);
  for (int i = 0; i < 10; ++i) {
    Eigen::Vector3d x = Eigen::Vector3d::Random();
    Eigen::Vector3d v1 = l1(x);
    Eigen::Vector3d v2 = l2(x);
    compDouble((lcross(x) - v1.cross(v2)).norm(), 0.);
  }
}

BOOST_AUTO_TEST_CASE(crossProductBezierLinearVariable) {
  bezier_linear_variable_t::t_point_t vec1;
  bezier_linear_variable_t::t_point_t vec2;
  bezier_linear_variable_t::t_point_t zeroVec;
  for (int i = 0; i < 3; ++i) {
    vec1.push_back(linear_variable_t(Eigen::Matrix3d::Identity() * i,
                                     Eigen::Vector3d::Random()));
    vec2.push_back(linear_variable_t(Eigen::Matrix3d::Identity() * i,
                                     Eigen::Vector3d::Random()));
    zeroVec.push_back(linear_variable_t());
  }
  for (int i = 0; i < 2; ++i) {
    vec1.push_back(linear_variable_t(Eigen::Matrix3d::Identity() * i,
                                     Eigen::Vector3d::Random()));
  }
  bezier_linear_variable_t p1(vec1.begin(), vec1.end(), 0., 1.);
  bezier_linear_variable_t p2(vec2.begin(), vec2.end(), 0., 1.);
  bezier_linear_variable_t pZero(zeroVec.begin(), zeroVec.end(), 0., 1.);
  bezier_linear_variable_t pCross(p1.cross(p2));
  bezier_linear_variable_t pCrossZero(p1.cross(pZero));
  bezier_linear_variable_t primitive =
      pCross.compute_primitive(1, *vec1.begin());
  for (double i = 0.; i <= 1.; ++i) {
    Eigen::Vector3d x = Eigen::Vector3d::Random();
    bezier_t fcross =
        evaluateLinear<bezier_t, bezier_linear_variable_t>(pCross, x);
    bezier_t fCrossZero =
        evaluateLinear<bezier_t, bezier_linear_variable_t>(pCrossZero, x);
    bezier_t f1 = evaluateLinear<bezier_t, bezier_linear_variable_t>(p1, x);
    bezier_t f2 = evaluateLinear<bezier_t, bezier_linear_variable_t>(p2, x);
    for (double i = 0.; i <= 10.; ++i) {
      double dt = i / 10.;
      Eigen::Vector3d v1 = f1(dt);
      Eigen::Vector3d v2 = f2(dt);
      compDouble((fcross(dt) - v1.cross(v2)).norm(), 0.);
      compDouble((fCrossZero(dt)).norm(), 0.);
    }
  }
}

BOOST_AUTO_TEST_CASE(polynomialPointOperations) {
  polynomial_t::coeff_t coeffs1 = Eigen::MatrixXd::Random(3, 5);
  Eigen::Vector3d point = Eigen::Vector3d::Random(3);
  polynomial_t p1(coeffs1, 0., 1.);

  polynomial_t pSum = p1 + point;
  polynomial_t pSumR = point + p1;
  polynomial_t pSub = p1 - point;
  polynomial_t pSubR = point - p1;
  polynomial_t pcross = p1.cross(point);
  for (double i = 0.; i <= 100.; ++i) {
    double dt = i / 100.;
    compDouble((pSum(dt) - (p1(dt) + point)).norm(), 0.);
    compDouble((pSumR(dt) - (p1(dt) + point)).norm(), 0.);
    compDouble((pSub(dt) - (p1(dt) - point)).norm(), 0.);
    compDouble((pSubR(dt) - (point - p1(dt))).norm(), 0.);
    Eigen::Vector3d p1dt = p1(dt);
    compDouble((pcross(dt) - (p1dt.cross(point))).norm(), 0.);
  }
}

BOOST_AUTO_TEST_CASE(polynomialOperations) {
  polynomial_t::coeff_t coeffs1 = Eigen::MatrixXd::Random(3, 5);
  polynomial_t::coeff_t coeffs2 = Eigen::MatrixXd::Random(3, 2);
  polynomial_t p1(coeffs1, 0., 1.);
  polynomial_t p2(coeffs2, 0., 1.);
  polynomial_t p3(coeffs2, 0., 0.5);
  polynomial_t p4(coeffs2, 0.1, 1.);
  polynomial_t p5(coeffs2, 0.1, .5);
  double k = 10.2;

  BOOST_CHECK_THROW(p1 + p3, std::exception);
  BOOST_CHECK_THROW(p1 - p3, std::exception);
  BOOST_CHECK_THROW(p1 + p4, std::exception);
  BOOST_CHECK_THROW(p1 - p4, std::exception);
  BOOST_CHECK_THROW(p1 + p5, std::exception);
  BOOST_CHECK_THROW(p1 - p5, std::exception);

  polynomial_t pSum = p1 + p2;
  polynomial_t pSumR = p2 + p1;
  polynomial_t pSub = p1 - p2;
  polynomial_t pSubR = p2 - p1;
  polynomial_t pdiv = p1 / k;
  polynomial_t pMul = p1 * k;
  polynomial_t pMulR = k * p1;
  polynomial_t pNeg = -p1;
  for (double i = 0.; i <= 100.; ++i) {
    double dt = i / 100.;
    compDouble((pSum(dt) - (p1(dt) + p2(dt))).norm(), 0.);
    compDouble((pSumR(dt) - (p1(dt) + p2(dt))).norm(), 0.);
    compDouble((pSub(dt) - (p1(dt) - p2(dt))).norm(), 0.);
    compDouble((pSubR(dt) - (p2(dt) - p1(dt))).norm(), 0.);
    compDouble((pMul(dt) - p1(dt) * k).norm(), 0.);
    compDouble((pMulR(dt) - p1(dt) * k).norm(), 0.);
    compDouble((pdiv(dt) - p1(dt) / k).norm(), 0.);
    compDouble((pNeg(dt) + p1(dt)).norm(), 0);
  }

  pSum = polynomial_t(p1);
  pSum += p2;
  pSub = p1;
  pSub -= p2;
  pdiv = p1;
  pdiv /= k;
  pMul = p1;
  pMul *= k;
  for (double i = 0.; i <= 100.; ++i) {
    double dt = i / 100.;
    compDouble((pSum(dt) - (p1(dt) + p2(dt))).norm(), 0.);
    compDouble((pSub(dt) - (p1(dt) - p2(dt))).norm(), 0.);
    compDouble((pMul(dt) - p1(dt) * k).norm(), 0.);
    compDouble((pdiv(dt) - p1(dt) / k).norm(), 0.);
  }
}

BOOST_AUTO_TEST_CASE(crossPoductPolynomials) {
  polynomial_t::coeff_t coeffs1 = Eigen::MatrixXd::Random(3, 5);
  polynomial_t::coeff_t coeffs2 = Eigen::MatrixXd::Random(3, 2);
  polynomial_t::coeff_t coeffsDim4 = Eigen::MatrixXd::Random(4, 2);
  polynomial_t p1(coeffs1, 0., 1.);
  polynomial_t p2(coeffs2, 0., 1.);
  polynomial_t p3(coeffs2, 0., 0.5);
  polynomial_t p4(coeffs2, 0.1, 1.);
  polynomial_t p5(coeffs2, 0.1, .5);
  polynomial_t pDim4(coeffsDim4, 0., 1.);

  // testing reduction
  BOOST_CHECK_EQUAL(p1.cross(p1).degree(), 0);

  BOOST_CHECK_THROW(p1.cross(p3), std::exception);
  BOOST_CHECK_THROW(p1.cross(p4), std::exception);
  BOOST_CHECK_THROW(p1.cross(p5), std::exception);
  BOOST_CHECK_THROW(p1.cross(pDim4), std::exception);
  BOOST_CHECK_THROW(pDim4.cross(p1), std::exception);

  polynomial_t pCross = p1.cross(p2);
  for (double i = 0.; i <= 100.; ++i) {
    double dt = i / 100.;
    Eigen::Vector3d v1 = p1(dt);
    Eigen::Vector3d v2 = p2(dt);
    compDouble((pCross(dt) - v1.cross(v2)).norm(), 0.);
  }
}

BOOST_AUTO_TEST_CASE(crossPoductPolynomialSimplification) {
  polynomial_t::coeff_t coeffs1 = Eigen::MatrixXd::Random(3, 5);
  polynomial_t::coeff_t coeffs2 = Eigen::MatrixXd::Random(3, 3);
  coeffs2.col(2) = coeffs1.col(4);
  polynomial_t p1(coeffs1, 0., 1.);
  polynomial_t p2(coeffs2, 0., 1.);

  polynomial_t pCross = p1.cross(p2);
  BOOST_CHECK_EQUAL(pCross.degree(), 5);
  for (double i = 0.; i <= 100.; ++i) {
    double dt = i / 100.;
    Eigen::Vector3d v1 = p1(dt);
    Eigen::Vector3d v2 = p2(dt);
    compDouble((pCross(dt) - v1.cross(v2)).norm(), 0.);
  }
}

BOOST_AUTO_TEST_SUITE_END()
