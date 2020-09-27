#define BOOST_TEST_MODULE test_operations

#include "curves/fwd.h"
#include "curves/exact_cubic.h"
#include "curves/bezier_curve.h"
#include "curves/polynomial.h"
#include "curves/helpers/effector_spline.h"
#include "curves/helpers/effector_spline_rotation.h"
#include "curves/curve_conversion.h"
#include "curves/cubic_hermite_spline.h"
#include "curves/piecewise_curve.h"
#include "curves/optimization/definitions.h"
#include "load_problem.h"
#include "curves/so3_linear.h"
#include "curves/se3_curve.h"
#include "curves/serialization/curves.hpp"
#include <boost/test/included/unit_test.hpp>

using namespace curves;

BOOST_AUTO_TEST_SUITE(BOOST_TEST_MODULE)

BOOST_AUTO_TEST_CASE(bezierOperations, * boost::unit_test::tolerance(0.001)) {
    t_pointX_t vec1;
    t_pointX_t vec2;
    for (int i =0; i<3; ++i)
    {
        vec1.push_back(Eigen::Vector3d::Random());
        vec2.push_back(Eigen::Vector3d::Random());
    }
    for (int i =0; i<2; ++i)
    {
        vec1.push_back(Eigen::Vector3d::Random());
    }
    bezier_t p1(vec1.begin(),vec1.end(),0.,1.);
    bezier_t p2(vec2.begin(),vec2.end(),0.,1.);
    bezier_t p3(vec2.begin(),vec2.end(),0.,0.5);
    bezier_t p4(vec2.begin(),vec2.end(),0.1,1.);
    bezier_t p5(vec2.begin(),vec2.end(),0.1,.5);
    double k = 10.2;

    BOOST_CHECK_THROW( p1 + p3, std::exception );
    BOOST_CHECK_THROW( p1 - p3, std::exception );
    BOOST_CHECK_THROW( p1 + p4, std::exception );
    BOOST_CHECK_THROW( p1 - p4, std::exception );
    BOOST_CHECK_THROW( p1 + p5, std::exception );
    BOOST_CHECK_THROW( p1 - p5, std::exception );

    bezier_t pSum  = p1 + p2;
    bezier_t pSumR = p2 + p1;
    bezier_t pSub  = p1 - p2;
    bezier_t pSubR = p2 - p1;
    bezier_t pdiv  = p1 / k;
    bezier_t pMul  = p1 * k;
    bezier_t pMulR = k  * p1;
    bezier_t pNeg  = -p1;
    for (double i = 0.; i <=100.; ++i ){
        double dt = i / 100.;
        BOOST_TEST(( pSum(dt) - (p1(dt)+p2(dt))).norm()==0.);
        BOOST_TEST((pSumR(dt) - (p1(dt)+p2(dt))).norm()==0.);
        BOOST_TEST(( pSub(dt) - (p1(dt)-p2(dt))).norm()==0.);
        BOOST_TEST((pSubR(dt) - (p2(dt)-p1(dt))).norm()==0.);
        BOOST_TEST(( pMul(dt) -  p1(dt)*k).norm()==0.);
        BOOST_TEST((pMulR(dt) -  p1(dt)*k).norm()==0.);
        BOOST_TEST(( pdiv(dt) -  p1(dt)/k).norm()==0.);
        BOOST_TEST(( pNeg(dt) +  p1(dt)).norm() == 0);
    }

    pSum = bezier_t(p1); pSum += p2;
    pSub = p1; pSub -= p2;
    pdiv = p1; pdiv /= k;
    pMul = p1; pMul *= k;
    for (double i = 0.; i <=100.; ++i ){
        double dt = i / 100.;
        BOOST_TEST(( pSum(dt) - (p1(dt)+p2(dt))).norm()==0.);
        BOOST_TEST(( pSub(dt) - (p1(dt)-p2(dt))).norm()==0.);
        BOOST_TEST(( pMul(dt) -  p1(dt)*k).norm()==0.);
        BOOST_TEST(( pdiv(dt) -  p1(dt)/k).norm()==0.);
    }
}

BOOST_AUTO_TEST_CASE(polynomialOperations, * boost::unit_test::tolerance(0.001)) {
    polynomial_t::coeff_t coeffs1 = Eigen::MatrixXd::Random(3,5);
    polynomial_t::coeff_t coeffs2 = Eigen::MatrixXd::Random(3,2);
    polynomial_t p1(coeffs1,0.,1.);
    polynomial_t p2(coeffs2,0.,1.);
    polynomial_t p3(coeffs2,0.,0.5);
    polynomial_t p4(coeffs2,0.1,1.);
    polynomial_t p5(coeffs2,0.1,.5);
    double k = 10.2;

    BOOST_CHECK_THROW( p1 + p3, std::exception );
    BOOST_CHECK_THROW( p1 - p3, std::exception );
    BOOST_CHECK_THROW( p1 + p4, std::exception );
    BOOST_CHECK_THROW( p1 - p4, std::exception );
    BOOST_CHECK_THROW( p1 + p5, std::exception );
    BOOST_CHECK_THROW( p1 - p5, std::exception );

    polynomial_t pSum  = p1 + p2;
    polynomial_t pSumR = p2 + p1;
    polynomial_t pSub  = p1 - p2;
    polynomial_t pSubR = p2 - p1;
    polynomial_t pdiv  = p1 / k;
    polynomial_t pMul  = p1 * k;
    polynomial_t pMulR = k  * p1;
    polynomial_t pNeg  = -p1;
    for (double i = 0.; i <=100.; ++i ){
        double dt = i / 100.;
        BOOST_TEST(( pSum(dt) - (p1(dt)+p2(dt))).norm()==0.);
        BOOST_TEST((pSumR(dt) - (p1(dt)+p2(dt))).norm()==0.);
        BOOST_TEST(( pSub(dt) - (p1(dt)-p2(dt))).norm()==0.);
        BOOST_TEST((pSubR(dt) - (p2(dt)-p1(dt))).norm()==0.);
        BOOST_TEST(( pMul(dt) -  p1(dt)*k).norm()==0.);
        BOOST_TEST((pMulR(dt) -  p1(dt)*k).norm()==0.);
        BOOST_TEST(( pdiv(dt) -  p1(dt)/k).norm()==0.);
        BOOST_TEST(( pNeg(dt) +  p1(dt)).norm() == 0);
    }

    pSum = polynomial_t(p1); pSum += p2;
    pSub = p1; pSub -= p2;
    pdiv = p1; pdiv /= k;
    pMul = p1; pMul *= k;
    for (double i = 0.; i <=100.; ++i ){
        double dt = i / 100.;
        BOOST_TEST(( pSum(dt) - (p1(dt)+p2(dt))).norm()==0.);
        BOOST_TEST(( pSub(dt) - (p1(dt)-p2(dt))).norm()==0.);
        BOOST_TEST(( pMul(dt) -  p1(dt)*k).norm()==0.);
        BOOST_TEST(( pdiv(dt) -  p1(dt)/k).norm()==0.);
    }
}

BOOST_AUTO_TEST_SUITE_END()
