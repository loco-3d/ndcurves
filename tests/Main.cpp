
#include "curves/exact_cubic.h"
#include "curves/bezier_curve.h"
#include "curves/polynomial.h"
#include "curves/helpers/effector_spline.h"
#include "curves/helpers/effector_spline_rotation.h"
#include "curves/curve_conversion.h"
#include "curves/cubic_hermite_spline.h"
#include "curves/piecewise_curve.h"

#include <string>
#include <iostream>
#include <cmath>

using namespace std;

namespace curves
{
typedef Eigen::Vector3d point_t;
typedef Eigen::Vector3d tangent_t;
typedef std::vector<point_t,Eigen::aligned_allocator<point_t> >  t_point_t;
typedef polynomial  <double, double, 3, true, point_t, t_point_t> polynomial_t;
typedef exact_cubic <double, double, 3, true, point_t> exact_cubic_t;
typedef bezier_curve  <double, double, 3, true, point_t> bezier_curve_t;
typedef exact_cubic_t::spline_constraints spline_constraints_t;
typedef std::pair<double, point_t> Waypoint;
typedef std::vector<Waypoint> T_Waypoint;


typedef Eigen::Matrix<double,1,1> point_one;
typedef polynomial<double, double, 1, true, point_one> polynom_one;
typedef exact_cubic   <double, double, 1, true, point_one> exact_cubic_one;
typedef std::pair<double, point_one> WaypointOne;
typedef std::vector<WaypointOne> T_WaypointOne;

typedef cubic_hermite_spline <double, double, 3, true, point_t> cubic_hermite_spline_t;
typedef std::pair<point_t, tangent_t> pair_point_tangent_t;
typedef std::vector<pair_point_tangent_t,Eigen::aligned_allocator<pair_point_tangent_t> > t_pair_point_tangent_t;

typedef piecewise_curve <double, double, 3, true, point_t, t_point_t, polynomial_t> piecewise_polynomial_curve_t;
typedef piecewise_curve <double, double, 3, true, point_t, t_point_t, bezier_curve_t> piecewise_bezier_curve_t;
typedef piecewise_curve <double, double, 3, true, point_t, t_point_t, cubic_hermite_spline_t> piecewise_cubic_hermite_curve_t;

const double margin = 1e-3;

bool QuasiEqual(const double a, const double b)
{
    return std::fabs(a-b)<margin;
}

} // namespace curves

using namespace curves;

ostream& operator<<(ostream& os, const point_t& pt)
{
    os << "(" << pt.x() << ", " << pt.y() << ", " << pt.z() << ")";
    return os;
}

void ComparePoints(const Eigen::VectorXd& pt1, const Eigen::VectorXd& pt2, const std::string& errmsg, bool& error, bool notequal = false)
{
    if(!QuasiEqual((pt1-pt2).norm(), 0.0)  && !notequal)
    {
        error = true;
        std::cout << errmsg << pt1.transpose() << " ; " << pt2.transpose() << std::endl;
    }
}

template<typename curve1, typename curve2>
void CompareCurves(curve1 c1, curve2 c2, const std::string& errMsg, bool& error)
{
    double T_min = c1.min();
    double T_max = c1.max();
    if (!QuasiEqual(T_min, c2.min()) || !QuasiEqual(T_max, c2.max()))
    {
        std::cout << "CompareCurves, ERROR, time min and max of curves do not match ["<<T_min<<","<<T_max<<"] " 
                << " and ["<<c2.min()<<","<<c2.max()<<"] "<<std::endl;
        error = true;
    }
    else
    {
        // derivative in T_min and T_max
        ComparePoints(c1.derivate(T_min,1),c2.derivate(T_min,1),errMsg, error, false);
        ComparePoints(c1.derivate(T_max,1),c2.derivate(T_max,1),errMsg, error, false);
        // Test values on curves
        for(double i =T_min; i<T_max; i+=0.02)
        {
            ComparePoints(c1(i),c2(i),errMsg, error, false);
            ComparePoints(c1(i),c2(i),errMsg, error, false);
        }
    }
}

/*Cubic Function tests*/
void CubicFunctionTest(bool& error)
{
    std::string errMsg("In test CubicFunctionTest ; unexpected result for x ");
    point_t a(1,2,3);
    point_t b(2,3,4);
    point_t c(3,4,5);
    point_t d(3,6,7);
    t_point_t vec;
    vec.push_back(a);
    vec.push_back(b);
    vec.push_back(c);
    vec.push_back(d);
    polynomial_t cf(vec.begin(), vec.end(), 0, 1);
    point_t res1;
    res1 =cf(0);
    point_t x0(1,2,3);
    ComparePoints(x0, res1, errMsg + "(0) ", error);

    point_t x1(9,15,19);
    res1 =cf(1);
    ComparePoints(x1, res1, errMsg + "(1) ", error);

    point_t x2(3.125,5.25,7.125);
    res1 =cf(0.5);
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
    try
    {
        cf2(0.4);
    }
    catch(...)
    {
        error = false;
    }
    if(error)
    {
        std::cout << "Evaluation of cubic cf2 error, 0.4 should be an out of range value\n";
    }
    error = true;
    try
    {
        cf2(1.1);
    }
    catch(...)
    {
        error = false;
    }
    if(error)
    {
        std::cout << "Evaluation of cubic cf2 error, 1.1 should be an out of range value\n";
    }
    if (!QuasiEqual(cf.max(), 1.0))
    {
        error = true;
        std::cout << "Evaluation of cubic cf error, MaxBound should be equal to 1\n";
    }
    if (!QuasiEqual(cf.min(), 0.0))
    {
        error = true;
        std::cout << "Evaluation of exactCubic error, MinBound should be equal to 0\n";
    }
}

/*bezier_curve Function tests*/
void BezierCurveTest(bool& error)
{
    std::string errMsg("In test BezierCurveTest ; unexpected result for x ");
    point_t a(1,2,3);
    point_t b(2,3,4);
    point_t c(3,4,5);
    point_t d(3,6,7);

    std::vector<point_t> params;
    params.push_back(a);

    // 1d curve in [0,1]
    bezier_curve_t cf1(params.begin(), params.end());
    point_t res1;
    res1 = cf1(0);
    point_t x10 = a ;
    ComparePoints(x10, res1, errMsg + "1(0) ", error);
    res1 =  cf1(1);
    ComparePoints(x10, res1, errMsg + "1(1) ", error);

    // 2d curve in [0,1]
    params.push_back(b);
    bezier_curve_t cf(params.begin(), params.end());
    res1 = cf(0);
    point_t x20 = a ;
    ComparePoints(x20, res1, errMsg + "2(0) ", error);

    point_t x21 = b;
    res1 = cf(1);
    ComparePoints(x21, res1, errMsg + "2(1) ", error);

    //3d curve in [0,1]
    params.push_back(c);
    bezier_curve_t cf3(params.begin(), params.end());
    res1 = cf3(0);
    ComparePoints(a, res1, errMsg + "3(0) ", error);

    res1 = cf3(1);
    ComparePoints(c, res1, errMsg + "3(1) ", error);

    //4d curve in [1,2]
    params.push_back(d);
    bezier_curve_t cf4(params.begin(), params.end(), 1., 2.);

    //testing bernstein polynomials
    bezier_curve_t cf5(params.begin(), params.end(),1.,2.);
    std::string errMsg2("In test BezierCurveTest ; Bernstein polynomials do not evaluate as analytical evaluation");
    for(double d = 1.; d <2.; d+=0.1)
    {
        ComparePoints( cf5.evalBernstein(d) , cf5 (d), errMsg2, error);
        ComparePoints( cf5.evalHorner(d) , cf5 (d), errMsg2, error);
        ComparePoints( cf5.compute_derivate(1).evalBernstein(d) , cf5.compute_derivate(1) (d), errMsg2, error);
        ComparePoints( cf5.compute_derivate(1).evalHorner(d) , cf5.compute_derivate(1) (d), errMsg2, error);
    }

    bool error_in(true);

    try
    {
        cf(-0.4);
    } catch(...)
    {
        error_in = false;
    }
    if(error_in)
    {
        std::cout << "Evaluation of bezier cf error, -0.4 should be an out of range value\n";
        error = true;
    }
    error_in = true;

    try
    {
        cf(1.1);
    } catch(...)
    {
        error_in = false;
    }
    if(error_in)
    {
        std::cout << "Evaluation of bezier cf error, 1.1 should be an out of range value\n";
        error = true;
    }
    if (!QuasiEqual(cf.max(),1.0))
    {
        error = true;
        std::cout << "Evaluation of bezier cf error, MaxBound should be equal to 1\n";
    }
    if (!QuasiEqual(cf.min(),0.0))
    {
        error = true;
        std::cout << "Evaluation of bezier cf error, MinBound should be equal to 1\n";
    }
}

#include <ctime>
void BezierCurveTestCompareHornerAndBernstein(bool&) // error
{
    using namespace std;
    std::vector<double> values;
    for (int i =0; i < 100000; ++i)
        values.push_back(rand()/RAND_MAX);

    //first compare regular evaluation (low dim pol)
    point_t a(1,2,3);
    point_t b(2,3,4);
    point_t c(3,4,5);
    point_t d(3,6,7);
    point_t e(3,61,7);
    point_t f(3,56,7);
    point_t g(3,36,7);
    point_t h(43,6,7);
    point_t i(3,6,77);

    std::vector<point_t> params;
    params.push_back(a);
    params.push_back(b);
    params.push_back(c);

    // 3d curve
    bezier_curve_t cf(params.begin(), params.end()); // defined in [0,1]

    // Check all evaluation of bezier curve
    clock_t s0,e0,s1,e1,s2,e2,s3,e3;
    s0 = clock();
    for(std::vector<double>::const_iterator cit = values.begin(); cit != values.end(); ++cit)
    {
        cf(*cit);
    }
    e0 = clock();

    s1 = clock();
    for(std::vector<double>::const_iterator cit = values.begin(); cit != values.end(); ++cit)
    {
        cf.evalBernstein(*cit);
    }
    e1 = clock();

    s2 = clock();
    for(std::vector<double>::const_iterator cit = values.begin(); cit != values.end(); ++cit)
    {
        cf.evalHorner(*cit);
    }
    e2 = clock();

    s3 = clock();
    for(std::vector<double>::const_iterator cit = values.begin(); cit != values.end(); ++cit)
    {
        cf.evalDeCasteljau(*cit);
    }
    e3 = clock();

    std::cout << "time for analytical eval  " <<   double(e0 - s0) / CLOCKS_PER_SEC << std::endl;
    std::cout << "time for bernstein eval   "  <<   double(e1 - s1) / CLOCKS_PER_SEC << std::endl;
    std::cout << "time for horner eval      "     <<   double(e2 - s2) / CLOCKS_PER_SEC << std::endl;
    std::cout << "time for deCasteljau eval "     <<   double(e3 - s3) / CLOCKS_PER_SEC << std::endl;

    std::cout << "now with high order polynomial "    << std::endl;

    params.push_back(d);
    params.push_back(e);
    params.push_back(f);
    params.push_back(g);
    params.push_back(h);
    params.push_back(i);

    bezier_curve_t cf2(params.begin(), params.end());

    s1 = clock();
    for(std::vector<double>::const_iterator cit = values.begin(); cit != values.end(); ++cit)
    {
        cf2.evalBernstein(*cit);
    }
    e1 = clock();

    s2 = clock();
    for(std::vector<double>::const_iterator cit = values.begin(); cit != values.end(); ++cit)
    {
        cf2.evalHorner(*cit);
    }
    e2 = clock();

    s0 = clock();
    for(std::vector<double>::const_iterator cit = values.begin(); cit != values.end(); ++cit)
    {
        cf2(*cit);
    }
    e0 = clock();

    s3 = clock();
    for(std::vector<double>::const_iterator cit = values.begin(); cit != values.end(); ++cit)
    {
        cf2.evalDeCasteljau(*cit);
    }
    e3 = clock();

    std::cout << "time for analytical eval  " <<   double(e0 - s0) / CLOCKS_PER_SEC << std::endl;
    std::cout << "time for bernstein eval   "  <<   double(e1 - s1) / CLOCKS_PER_SEC << std::endl;
    std::cout << "time for horner eval      "     <<   double(e2 - s2) / CLOCKS_PER_SEC << std::endl;
    std::cout << "time for deCasteljau eval "     <<   double(e3 - s3) / CLOCKS_PER_SEC << std::endl;

}

void BezierDerivativeCurveTest(bool& error)
{
    std::string errMsg("In test BezierDerivativeCurveTest ;, Error While checking value of point on curve : ");
    point_t a(1,2,3);
    point_t b(2,3,4);
    point_t c(3,4,5);

    std::vector<point_t> params;
    params.push_back(a);
    params.push_back(b);

    params.push_back(c);
    bezier_curve_t cf3(params.begin(), params.end());

    ComparePoints(cf3(0), cf3.derivate(0.,0), errMsg, error);
    ComparePoints(cf3(0), cf3.derivate(0.,1), errMsg, error, true);
    ComparePoints(point_t::Zero(), cf3.derivate(0.,100), errMsg, error);
}

void BezierDerivativeCurveTimeReparametrizationTest(bool& error)
{
    std::string errMsg("In test BezierDerivativeCurveTimeReparametrizationTest, Error While checking value of point on curve : ");
    point_t a(1,2,3);
    point_t b(2,3,4);
    point_t c(3,4,5);
    point_t d(3,4,5);
    point_t e(3,4,5);
    point_t f(3,4,5);

    std::vector<point_t> params;
    params.push_back(a);
    params.push_back(b);
    params.push_back(c);
    params.push_back(d);
    params.push_back(e);
    params.push_back(f);
    double Tmin = 0.;
    double Tmax = 2.;
    double diffT = Tmax-Tmin;
    bezier_curve_t cf(params.begin(), params.end());
    bezier_curve_t cfT(params.begin(), params.end(),Tmin,Tmax);

    ComparePoints(cf(0.5), cfT(1), errMsg, error);
    ComparePoints(cf.derivate(0.5,1), cfT.derivate(1,1) * (diffT), errMsg, error);
    ComparePoints(cf.derivate(0.5,2), cfT.derivate(1,2) * diffT*diffT, errMsg, error);
}

void BezierDerivativeCurveConstraintTest(bool& error)
{
    std::string errMsg0("In test BezierDerivativeCurveConstraintTest, Error While checking value of point on curve : ");
    point_t a(1,2,3);
    point_t b(2,3,4);
    point_t c(3,4,5);

    bezier_curve_t::curve_constraints_t constraints;
    constraints.init_vel = point_t(-1,-1,-1);
    constraints.init_acc = point_t(-2,-2,-2);
    constraints.end_vel = point_t(-10,-10,-10);
    constraints.end_acc = point_t(-20,-20,-20);

    std::vector<point_t> params;
    params.push_back(a);
    params.push_back(b);
    params.push_back(c);

    bezier_curve_t::num_t T_min = 1.0;
    bezier_curve_t::num_t T_max = 3.0;
    bezier_curve_t cf(params.begin(), params.end(), constraints, T_min, T_max);

    ComparePoints(a, cf(T_min), errMsg0, error);
    ComparePoints(c, cf(T_max), errMsg0, error);
    ComparePoints(constraints.init_vel, cf.derivate(T_min,1), errMsg0, error);
    ComparePoints(constraints.end_vel , cf.derivate(T_max,1), errMsg0, error);
    ComparePoints(constraints.init_acc, cf.derivate(T_min,2), errMsg0, error);
    ComparePoints(constraints.end_vel, cf.derivate(T_max,1), errMsg0, error);
    ComparePoints(constraints.end_acc, cf.derivate(T_max,2), errMsg0, error);

    std::string errMsg1("In test BezierDerivativeCurveConstraintTest, Error While checking checking degree of bezier curve :");
    std::string errMsg2("In test BezierDerivativeCurveConstraintTest, Error While checking checking size of bezier curve :");
    if (cf.degree_ != params.size() + 3)
    {
        error = true;
        std::cout << errMsg1 << cf.degree_ << " ; " << params.size()+3 << std::endl;
    }
    if (cf.size_   != params.size() + 4)
    {
        error = true;
        std::cout << errMsg2 << cf.size_ << " ; " << params.size()+4 << std::endl;
    }
}


void toPolynomialConversionTest(bool& error)
{
    // bezier to polynomial
    std::string errMsg("In test BezierToPolynomialConversionTest, Error While checking value of point on curve : ");
    point_t a(1,2,3);
    point_t b(2,3,4);
    point_t c(3,4,5);
    point_t d(3,6,7);
    point_t e(3,61,7);
    point_t f(3,56,7);
    point_t g(3,36,7);
    point_t h(43,6,7);
    point_t i(3,6,77);

    std::vector<point_t> control_points;
    control_points.push_back(a);
    control_points.push_back(b);
    control_points.push_back(c);
    control_points.push_back(d);
    control_points.push_back(e);
    control_points.push_back(f);
    control_points.push_back(g);
    control_points.push_back(h);
    control_points.push_back(i);

    bezier_curve_t::num_t T_min = 1.0;
    bezier_curve_t::num_t T_max = 3.0;
    bezier_curve_t bc(control_points.begin(), control_points.end(),T_min, T_max);
    polynomial_t pol = polynomial_from_curve<polynomial_t, bezier_curve_t>(bc);
    CompareCurves<polynomial_t, bezier_curve_t>(pol, bc, errMsg, error);
}

void cubicConversionTest(bool& error)
{
    std::string errMsg0("In test CubicConversionTest - convert hermite to, Error While checking value of point on curve : ");
    std::string errMsg1("In test CubicConversionTest - convert bezier to, Error While checking value of point on curve : ");
    std::string errMsg2("In test CubicConversionTest - convert polynomial to, Error While checking value of point on curve : ");
    
    // Create cubic hermite spline : Test hermite to bezier/polynomial
    point_t p0(1,2,3);
    point_t m0(2,3,4);
    point_t p1(3,4,5);
    point_t m1(3,6,7);
    pair_point_tangent_t pair0(p0,m0);
    pair_point_tangent_t pair1(p1,m1);
    t_pair_point_tangent_t control_points;
    control_points.push_back(pair0);
    control_points.push_back(pair1);
    std::vector< double > time_control_points;
    polynomial_t::num_t T_min = 1.0;
    polynomial_t::num_t T_max = 3.0;
    time_control_points.push_back(T_min);
    time_control_points.push_back(T_max);
    cubic_hermite_spline_t chs0(control_points.begin(), control_points.end(), time_control_points);
    
    // hermite to bezier
    //std::cout<<"======================= \n";
    //std::cout<<"hermite to bezier \n";
    bezier_curve_t bc0 = bezier_from_curve<bezier_curve_t, cubic_hermite_spline_t>(chs0);
    CompareCurves<cubic_hermite_spline_t, bezier_curve_t>(chs0, bc0, errMsg0, error);

    // hermite to pol
    //std::cout<<"======================= \n";
    //std::cout<<"hermite to polynomial \n";
    polynomial_t pol0 = polynomial_from_curve<polynomial_t, cubic_hermite_spline_t>(chs0);
    CompareCurves<cubic_hermite_spline_t, polynomial_t>(chs0, pol0, errMsg0, error);

    // pol to hermite
    //std::cout<<"======================= \n";
    //std::cout<<"polynomial to hermite \n";
    cubic_hermite_spline_t chs1 = hermite_from_curve<cubic_hermite_spline_t, polynomial_t>(pol0);
    CompareCurves<polynomial_t, cubic_hermite_spline_t>(pol0,chs1,errMsg2,error);

    // pol to bezier
    //std::cout<<"======================= \n";
    //std::cout<<"polynomial to bezier \n";
    bezier_curve_t bc1 = bezier_from_curve<bezier_curve_t, polynomial_t>(pol0);
    CompareCurves<bezier_curve_t, polynomial_t>(bc1, pol0, errMsg2, error);

    // Bezier to pol
    //std::cout<<"======================= \n";
    //std::cout<<"bezier to polynomial \n";
    polynomial_t pol1 = polynomial_from_curve<polynomial_t, bezier_curve_t>(bc0);
    CompareCurves<bezier_curve_t, polynomial_t>(bc0, pol1, errMsg1, error);

    // bezier => hermite
    //std::cout<<"======================= \n";
    //std::cout<<"bezier to hermite \n";
    cubic_hermite_spline_t chs2 = hermite_from_curve<cubic_hermite_spline_t, bezier_curve_t>(bc0);
    CompareCurves<bezier_curve_t, cubic_hermite_spline_t>(bc0, chs2, errMsg1, error);
}

/*Exact Cubic Function tests*/
void ExactCubicNoErrorTest(bool& error)
{
    // Create an exact cubic spline with 7 waypoints => 6 polynomials defined in [0.0,3.0]
    curves::T_Waypoint waypoints;
    for(double i = 0.0; i <= 3.0; i = i + 0.5)
    {
        waypoints.push_back(std::make_pair(i,point_t(i,i,i)));
    }
    exact_cubic_t exactCubic(waypoints.begin(), waypoints.end());

    // Test number of polynomials in exact cubic
    std::size_t numberSegments = exactCubic.getNumberSplines();
    if (numberSegments != 6)
    {
        error = true;
        std::cout << "In ExactCubicNoErrorTest, Error While checking number of splines" << 
                    numberSegments << " ; " << 6 << std::endl;
    }

    // Test getSplineAt function
    for (std::size_t i=0; i<numberSegments; i++)
    {
        exactCubic.getSplineAt(i);
    }

    // Other tests
    try
    {
        exactCubic(0.0);
        exactCubic(3.0);
    }
    catch(...)
    {
        error = true;
        std::cout << "Evaluation of ExactCubicNoErrorTest error when testing value on bounds\n";
    }
    error = true;
    try
    {
        exactCubic(3.2);
    }
    catch(...)
    {
        error = false;
    }
    if(error)
    {
        std::cout << "Evaluation of exactCubic cf error, 3.2 should be an out of range value\n";
    }

    if (!QuasiEqual(exactCubic.max(),3.0))
    {
        error = true;
        std::cout << "Evaluation of exactCubic error, MaxBound should be equal to 3\n";
    }
    if (!QuasiEqual(exactCubic.min(),0.0))
    {
        error = true;
        std::cout << "Evaluation of exactCubic error, MinBound should be equal to 0\n";
    }
}

/*Exact Cubic Function tests*/
void ExactCubicTwoPointsTest(bool& error)
{
    // Create an exact cubic spline with 2 waypoints => 1 polynomial defined in [0.0,1.0]
    curves::T_Waypoint waypoints;
    for(double i = 0.0; i < 2.0; ++i)
    {
        waypoints.push_back(std::make_pair(i,point_t(i,i,i)));
    }
    exact_cubic_t exactCubic(waypoints.begin(), waypoints.end());

    point_t res1 = exactCubic(0);
    std::string errmsg0("in ExactCubicTwoPointsTest, Error While checking that given wayPoints  are crossed (expected / obtained)");
    ComparePoints(point_t(0,0,0), res1, errmsg0, error);

    res1 = exactCubic(1);
    ComparePoints(point_t(1,1,1), res1, errmsg0, error);

    // Test number of polynomials in exact cubic
    std::size_t numberSegments = exactCubic.getNumberSplines();
    if (numberSegments != 1)
    {
        error = true;
        std::cout << "In ExactCubicTwoPointsTest, Error While checking number of splines" << 
                    numberSegments << " ; " << 1 << std::endl;
    }

    // Test getSplineAt
    std::string errmsg1("in ExactCubicTwoPointsTest, Error While checking value on curve");
    ComparePoints(exactCubic(0.5), (exactCubic.getSplineAt(0))(0.5), errmsg1, error);
}

void ExactCubicOneDimTest(bool& error)
{
    curves::T_WaypointOne waypoints;
    point_one zero; zero(0,0) = 9;
    point_one one; one(0,0) = 14;
    point_one two; two(0,0) = 25;
    waypoints.push_back(std::make_pair(0., zero));
    waypoints.push_back(std::make_pair(1., one));
    waypoints.push_back(std::make_pair(2., two));
    exact_cubic_one exactCubic(waypoints.begin(), waypoints.end());

    point_one res1 = exactCubic(0);
    std::string errmsg("in ExactCubicOneDim Error While checking that given wayPoints  are crossed (expected / obtained)");
    ComparePoints(zero, res1, errmsg, error);

    res1 = exactCubic(1);
    ComparePoints(one, res1, errmsg, error);
}

void CheckWayPointConstraint(const std::string& errmsg, const double step, const curves::T_Waypoint&, const exact_cubic_t* curve, bool& error )
{
    point_t res1;
    for(double i = 0; i <= 1; i = i + step)
    {
        res1 = (*curve)(i);
        ComparePoints(point_t(i,i,i), res1, errmsg, error);
    }
}

void CheckDerivative(const std::string& errmsg, const double eval_point, const std::size_t order, const point_t& target, const exact_cubic_t* curve, bool& error )
{
    point_t res1 = curve->derivate(eval_point,order);
    ComparePoints(target, res1, errmsg, error);
}


void ExactCubicPointsCrossedTest(bool& error)
{
    curves::T_Waypoint waypoints;
    for(double i = 0; i <= 1; i = i + 0.2)
    {
        waypoints.push_back(std::make_pair(i,point_t(i,i,i)));
    }
    exact_cubic_t exactCubic(waypoints.begin(), waypoints.end());
    std::string errmsg("Error While checking that given wayPoints are crossed (expected / obtained)");
    CheckWayPointConstraint(errmsg, 0.2, waypoints, &exactCubic, error);

}

void ExactCubicVelocityConstraintsTest(bool& error)
{
    curves::T_Waypoint waypoints;
    for(double i = 0; i <= 1; i = i + 0.2)
    {
        waypoints.push_back(std::make_pair(i,point_t(i,i,i)));
    }
    std::string errmsg("Error in ExactCubicVelocityConstraintsTest (1); while checking that given wayPoints are crossed (expected / obtained)");
    spline_constraints_t constraints;
    constraints.end_vel = point_t(0,0,0);
    constraints.init_vel = point_t(0,0,0);
    constraints.end_acc = point_t(0,0,0);
    constraints.init_acc = point_t(0,0,0);
    exact_cubic_t exactCubic(waypoints.begin(), waypoints.end(), constraints);
    // now check that init and end velocity are 0
    CheckWayPointConstraint(errmsg, 0.2, waypoints, &exactCubic, error);
    std::string errmsg3("Error in ExactCubicVelocityConstraintsTest (2); while checking derivative (expected / obtained)");
    // now check derivatives
    CheckDerivative(errmsg3,0,1,constraints.init_vel,&exactCubic, error);
    CheckDerivative(errmsg3,1,1,constraints.end_vel,&exactCubic, error);
    CheckDerivative(errmsg3,0,2,constraints.init_acc,&exactCubic, error);
    CheckDerivative(errmsg3,1,2,constraints.end_acc,&exactCubic, error);

    constraints.end_vel = point_t(1,2,3);
    constraints.init_vel = point_t(-1,-2,-3);
    constraints.end_acc = point_t(4,5,6);
    constraints.init_acc = point_t(-4,-4,-6);
    std::string errmsg2("Error in ExactCubicVelocityConstraintsTest (3); while checking that given wayPoints are crossed (expected / obtained)");
    exact_cubic_t exactCubic2(waypoints.begin(), waypoints.end(),constraints);
    CheckWayPointConstraint(errmsg2, 0.2, waypoints, &exactCubic2, error);

    std::string errmsg4("Error in ExactCubicVelocityConstraintsTest (4); while checking derivative (expected / obtained)");
    // now check derivatives
    CheckDerivative(errmsg4,0,1,constraints.init_vel,&exactCubic2, error);
    CheckDerivative(errmsg4,1,1,constraints.end_vel ,&exactCubic2, error);
    CheckDerivative(errmsg4,0,2,constraints.init_acc,&exactCubic2, error);
    CheckDerivative(errmsg4,1,2,constraints.end_acc ,&exactCubic2, error);
}

void CheckPointOnline(const std::string& errmsg, const point_t& A, const point_t& B, const double target, const exact_cubic_t* curve, bool& error )
{
    point_t res1 = curve->operator ()(target);
    point_t ar =(res1-A); ar.normalize();
    point_t rb =(B-res1); rb.normalize();
    if(ar.dot(rb) < 0.99999)
    {
        error = true;
        std::cout << errmsg << " ; " << A.transpose() << "\n ; " << B.transpose() << "\n ; " <<
                     target << " ; " << res1.transpose() <<  std::endl;
    }
}

void EffectorTrajectoryTest(bool& error)
{
    // create arbitrary trajectory
    curves::T_Waypoint waypoints;
    for(double i = 0; i <= 10; i = i + 2)
    {
        waypoints.push_back(std::make_pair(i,point_t(i,i,i)));
    }
    helpers::exact_cubic_t* eff_traj = helpers::effector_spline(waypoints.begin(),waypoints.end(),
                                                               Eigen::Vector3d::UnitZ(),Eigen::Vector3d(0,0,2),
                                                               1,0.02,1,0.5);
    point_t zero(0,0,0);
    point_t off1(0,0,1);
    point_t off2(10,10,10.02);
    point_t end(10,10,10);
    std::string errmsg("Error in EffectorTrajectoryTest; while checking waypoints (expected / obtained)");
    std::string errmsg2("Error in EffectorTrajectoryTest; while checking derivative (expected / obtained)");
    //first check start / goal positions
    ComparePoints(zero, (*eff_traj)(0), errmsg, error);
    ComparePoints(off1, (*eff_traj)(1), errmsg, error);
    ComparePoints(off2, (*eff_traj)(9.5), errmsg, error);
    ComparePoints(end , (*eff_traj)(10), errmsg, error);

    //then check offset at start / goal positions
    // now check derivatives
    CheckDerivative(errmsg2,0,1,zero,eff_traj, error);
    CheckDerivative(errmsg2,10,1,zero ,eff_traj, error);
    CheckDerivative(errmsg2,0,2,zero,eff_traj, error);
    CheckDerivative(errmsg2,10,2,zero ,eff_traj, error);

    //check that end and init splines are line
    std::string errmsg3("Error in EffectorTrajectoryTest; while checking that init/end splines are line (point A/ point B, time value / point obtained) \n");
    for(double i = 0.1; i<1; i+=0.1)
    {
        CheckPointOnline(errmsg3,(*eff_traj)(0),(*eff_traj)(1),i,eff_traj,error);
    }

    for(double i = 9.981; i<10; i+=0.002)
    {
        CheckPointOnline(errmsg3,(*eff_traj)(9.5),(*eff_traj)(10),i,eff_traj,error);
    }
    delete eff_traj;
}

helpers::quat_t GetXRotQuat(const double theta)
{
    Eigen::AngleAxisd m (theta, Eigen::Vector3d::UnitX());
    return helpers::quat_t(Eigen::Quaterniond(m).coeffs().data());
}

double GetXRotFromQuat(helpers::quat_ref_const_t q)
{
    Eigen::Quaterniond quat (q.data());
    Eigen::AngleAxisd m (quat);
    return m.angle() / M_PI * 180.;
}

void EffectorSplineRotationNoRotationTest(bool& error)
{
    // create arbitrary trajectory
    curves::T_Waypoint waypoints;
    for(double i = 0; i <= 10; i = i + 2)
    {
        waypoints.push_back(std::make_pair(i,point_t(i,i,i)));
    }
    helpers::effector_spline_rotation eff_traj(waypoints.begin(),waypoints.end());
    helpers::config_t q_init; q_init    << 0.,0.,0.,0.,0.,0.,1.;
    helpers::config_t q_end; q_end      << 10.,10.,10.,0.,0.,0.,1.;
    helpers::config_t q_to; q_to        << 0.,0,0.02,0.,0.,0.,1.;
    helpers::config_t q_land; q_land    << 10,10, 10.02, 0, 0.,0.,1.;
    helpers::config_t q_mod; q_mod      << 6.,6.,6.,0.,0.,0.,1.;
    std::string errmsg("Error in EffectorSplineRotationNoRotationTest; while checking waypoints (expected / obtained)");
    ComparePoints(q_init , eff_traj(0),    errmsg,error);
    ComparePoints(q_to   , eff_traj(0.02), errmsg,error);
    ComparePoints(q_land , eff_traj(9.98), errmsg,error);
    ComparePoints(q_mod  , eff_traj(6),    errmsg,error);
    ComparePoints(q_end  , eff_traj(10),   errmsg,error);
}

void EffectorSplineRotationRotationTest(bool& error)
{
    // create arbitrary trajectory
    curves::T_Waypoint waypoints;
    for(double i = 0; i <= 10; i = i + 2)
    {
        waypoints.push_back(std::make_pair(i,point_t(i,i,i)));
    }
    helpers::quat_t init_quat = GetXRotQuat(M_PI);
    helpers::effector_spline_rotation eff_traj(waypoints.begin(),waypoints.end(), init_quat);
    helpers::config_t q_init =  helpers::config_t::Zero(); q_init.tail<4>() = init_quat;
    helpers::config_t q_end; q_end      << 10.,10.,10.,0.,0.,0.,1.;
    helpers::config_t q_to   = q_init; q_to(2)  +=0.02;
    helpers::config_t q_land = q_end ; q_land(2)+=0.02;
    helpers::quat_t q_mod = GetXRotQuat(M_PI_2);;
    std::string errmsg("Error in EffectorSplineRotationRotationTest; while checking waypoints (expected / obtained)");
    ComparePoints(q_init, eff_traj(0),           errmsg,error);
    ComparePoints(q_to  , eff_traj(0.02),        errmsg,error);
    ComparePoints(q_land, eff_traj(9.98),        errmsg,error);
    ComparePoints(q_mod , eff_traj(5).tail<4>(), errmsg,error);
    ComparePoints(q_end , eff_traj(10),          errmsg,error);
}

void EffectorSplineRotationWayPointRotationTest(bool& error)
{
    // create arbitrary trajectory
    curves::T_Waypoint waypoints;
    for(double i = 0; i <= 10; i = i + 2)
    {
        waypoints.push_back(std::make_pair(i,point_t(i,i,i)));
    }
    helpers::quat_t init_quat = GetXRotQuat(0);
    helpers::t_waypoint_quat_t quat_waypoints_;


    helpers::quat_t q_pi_0 = GetXRotQuat(0);
    helpers::quat_t q_pi_2 = GetXRotQuat(M_PI_2);
    helpers::quat_t q_pi   = GetXRotQuat(M_PI);

    quat_waypoints_.push_back(std::make_pair(0.4,q_pi_0));
    quat_waypoints_.push_back(std::make_pair(6,q_pi_2));
    quat_waypoints_.push_back(std::make_pair(8,q_pi));


    helpers::effector_spline_rotation eff_traj(waypoints.begin(),waypoints.end(),
                                               quat_waypoints_.begin(), quat_waypoints_.end());
    helpers::config_t q_init =  helpers::config_t::Zero(); q_init.tail<4>() = init_quat;
    helpers::config_t q_end; q_end      << 10.,10.,10.,0.,0.,0.,1.; q_end.tail<4>() = q_pi;
    helpers::config_t q_mod; q_mod.head<3>() = point_t(6,6,6) ; q_mod.tail<4>() = q_pi_2;
    helpers::config_t q_to   = q_init; q_to(2)  +=0.02;
    helpers::config_t q_land = q_end ; q_land(2)+=0.02;
    std::string errmsg("Error in EffectorSplineRotationWayPointRotationTest; while checking waypoints (expected / obtained)");
    ComparePoints(q_init, eff_traj(0),           errmsg,error);
    ComparePoints(q_to  , eff_traj(0.02),        errmsg,error);
    ComparePoints(q_land, eff_traj(9.98),        errmsg,error);
    ComparePoints(q_mod , eff_traj(6), errmsg,error);
    ComparePoints(q_end , eff_traj(10),          errmsg,error);
}

void TestReparametrization(bool& error)
{
    helpers::rotation_spline s;
    const helpers::exact_cubic_constraint_one_dim& sp = s.time_reparam_;
    if (!QuasiEqual(sp.min(),0.0))
    {
        std::cout << "in TestReparametrization; min value is not 0, got " << sp.min() << std::endl;
        error = true;
    }
    if (!QuasiEqual(sp.max(),1.0))
    {
        std::cout << "in TestReparametrization; max value is not 1, got " << sp.max() << std::endl;
        error = true;
    }
    if (!QuasiEqual(sp(1)[0],1.0))
    {
        std::cout << "in TestReparametrization; end value is not 1, got " << sp(1)[0] << std::endl;
        error = true;
    }
    if (!QuasiEqual(sp(0)[0],0.0))
    {
        std::cout << "in TestReparametrization; init value is not 0, got " << sp(0)[0] << std::endl;
        error = true;
    }
    for(double i =0; i<1; i+=0.002)
    {
        if(sp(i)[0]>sp(i+0.002)[0])
        {
            std::cout << "in TestReparametrization; reparametrization not monotonous " << sp.max() << std::endl;
            error = true;
        }
    }
}

point_t randomPoint(const double min, const double max )
{
    point_t p;
    for(size_t i = 0 ; i < 3 ; ++i)
    {
        p[i] =  (rand()/(double)RAND_MAX ) * (max-min) + min;
    }
    return p;
}

void BezierEvalDeCasteljau(bool& error)
{
    using namespace std;
    std::vector<double> values;
    for (int i =0; i < 100000; ++i)
    {
        values.push_back(rand()/RAND_MAX);
    }

    //first compare regular evaluation (low dim pol)
    point_t a(1,2,3);
    point_t b(2,3,4);
    point_t c(3,4,5);
    point_t d(3,6,7);
    point_t e(3,61,7);
    point_t f(3,56,7);
    point_t g(3,36,7);
    point_t h(43,6,7);
    point_t i(3,6,77);

    std::vector<point_t> params;
    params.push_back(a);
    params.push_back(b);
    params.push_back(c);

    // 3d curve
    bezier_curve_t cf(params.begin(), params.end());

    std::string errmsg("Error in BezierEvalDeCasteljau; while comparing actual bezier evaluation and de Casteljau : ");
    for(std::vector<double>::const_iterator cit = values.begin(); cit != values.end(); ++cit)
    {
        ComparePoints(cf.evalDeCasteljau(*cit), cf(*cit), errmsg, error);
    }

    params.push_back(d);
    params.push_back(e);
    params.push_back(f);
    params.push_back(g);
    params.push_back(h);
    params.push_back(i);

    bezier_curve_t cf2(params.begin(), params.end());

    for(std::vector<double>::const_iterator cit = values.begin(); cit != values.end(); ++cit)
    {
        ComparePoints(cf.evalDeCasteljau(*cit), cf(*cit), errmsg, error);
    }

}

/**
 * @brief BezierSplitCurve test the 'split' method of bezier curve
 * @param error
 */
void BezierSplitCurve(bool& error)
{
    // test for degree 5
    size_t n = 5;
    double t_min = 0.2;
    double t_max = 10;
    double aux0, aux1;
    for(size_t i = 0 ; i < 1 ; ++i)
    {
        // build a random curve and split it at random time :
        //std::cout<<"build a random curve"<<std::endl;
        point_t a;
        std::vector<point_t> wps;
        for(size_t j = 0 ; j <= n ; ++j)
        {
            wps.push_back(randomPoint(-10.,10.));
        }
        double t0 = (rand()/(double)RAND_MAX )*(t_max-t_min) + t_min;
        double t1 = (rand()/(double)RAND_MAX )*(t_max-t0) + t0;
        double ts = (rand()/(double)RAND_MAX )*(t1-t0)+t0;

        bezier_curve_t c(wps.begin(), wps.end(),t0, t1);
        std::pair<bezier_curve_t,bezier_curve_t> cs = c.split(ts);

        // test on splitted curves :
        if(! ((c.degree_ == cs.first.degree_) && (c.degree_ == cs.second.degree_) ))
        {
            error = true;
            std::cout<<"BezierSplitCurve, ERROR Degree of the splitted curve are not the same as the original curve"<<std::endl;
        }
        aux0 = c.max()-c.min();
        aux1 = (cs.first.max()-cs.first.min() + cs.second.max()-cs.second.min());
        if(!QuasiEqual(aux0, aux1))
        {
            error = true;
            std::cout<<"BezierSplitCurve, ERROR duration of the splitted curve doesn't correspond to the original"<<std::endl;
        }
        if(!QuasiEqual(cs.first.max(), ts))
        {
            error = true;
            std::cout<<"BezierSplitCurve, ERROR timing of the splitted curve doesn't correspond to the original"<<std::endl;
        }

        std::string errmsg("BezierSplitCurve, ERROR initial point of the splitted curve doesn't correspond to the original");
        ComparePoints(c(t0), cs.first(t0), errmsg, error);
        errmsg = "BezierSplitCurve, ERROR splitting point of the splitted curve doesn't correspond to the original";
        ComparePoints(cs.first(ts), cs.second(ts), errmsg, error);
        errmsg = "BezierSplitCurve, ERROR final point of the splitted curve doesn't correspond to the original";
        ComparePoints(c(t1), cs.second(cs.second.max()), errmsg, error);

        // check along curve :
        double ti = t0;
        errmsg = "BezierSplitCurve, ERROR while checking value on curve and curves splitted";
        while(ti <= ts)
        {
            ComparePoints(cs.first(ti), c(ti), errmsg, error);
            ti += 0.01;
        }
        while(ti <= t1)
        {
            ComparePoints(cs.second(ti), c(ti), errmsg, error);
            ti += 0.01;
        }
    }
}


/* cubic hermite spline function test */
void CubicHermitePairsPositionDerivativeTest(bool& error)
{
    std::string errmsg1("in Cubic Hermite 2 pairs (pos,vel), Error While checking that given wayPoints are crossed (expected / obtained) : ");
    std::string errmsg2("in Cubic Hermite 2 points, Error While checking value of point on curve : ");
    std::string errmsg3("in Cubic Hermite 2 points, Error While checking value of tangent on curve : ");
    std::vector< pair_point_tangent_t > control_points;
    point_t res1;

    point_t p0(0.,0.,0.);
    point_t p1(1.,2.,3.);
    point_t p2(4.,4.,4.);

    tangent_t t0(0.5,0.5,0.5);
    tangent_t t1(0.1,0.2,-0.5);
    tangent_t t2(0.1,0.2,0.3);

    std::vector< double > time_control_points;

    // Two pairs
    control_points.clear();
    control_points.push_back(pair_point_tangent_t(p0,t0));
    control_points.push_back(pair_point_tangent_t(p1,t1));
    time_control_points.push_back(1.);  // Time at P0
    time_control_points.push_back(3.);  // Time at P1
    // Create cubic hermite spline
    cubic_hermite_spline_t cubic_hermite_spline_1Pair(control_points.begin(), control_points.end(), time_control_points);
    //Check
    res1 = cubic_hermite_spline_1Pair(1.);   // t=0
    ComparePoints(p0, res1, errmsg1, error);
    res1 = cubic_hermite_spline_1Pair(3.);   // t=1
    ComparePoints(p1, res1, errmsg1, error);
    // Test derivative : two pairs
    res1 = cubic_hermite_spline_1Pair.derivate(1.,1);
    ComparePoints(t0, res1, errmsg3, error);
    res1 = cubic_hermite_spline_1Pair.derivate(3.,1);
    ComparePoints(t1, res1, errmsg3, error);

    // Three pairs
    control_points.push_back(pair_point_tangent_t(p2,t2));
    time_control_points.clear();
    time_control_points.push_back(0.);  // Time at P0
    time_control_points.push_back(2.);  // Time at P1
    time_control_points.push_back(5.);  // Time at P2
    cubic_hermite_spline_t cubic_hermite_spline_2Pairs(control_points.begin(), control_points.end(), time_control_points);
    //Check
    res1 = cubic_hermite_spline_2Pairs(0.);  // t=0
    ComparePoints(p0, res1, errmsg1, error);
    res1 = cubic_hermite_spline_2Pairs(2.);  // t=2
    ComparePoints(p1, res1, errmsg2, error);
    res1 = cubic_hermite_spline_2Pairs(5.);  // t=5
    ComparePoints(p2, res1, errmsg1, error);
    // Test derivative : three pairs
    res1 = cubic_hermite_spline_2Pairs.derivate(0.,1);
    ComparePoints(t0, res1, errmsg3, error);
    res1 = cubic_hermite_spline_2Pairs.derivate(2.,1);
    ComparePoints(t1, res1, errmsg3, error);
    res1 = cubic_hermite_spline_2Pairs.derivate(5.,1);
    ComparePoints(t2, res1, errmsg3, error);
    // Test time control points by default [0,1] => with N control points : 
    // Time at P0= 0. | Time at P1= 1.0/(N-1) | Time at P2= 2.0/(N-1) | ... | Time at P_(N-1)= (N-1)/(N-1)= 1.0
    time_control_points.clear();
    time_control_points.push_back(0.);  // Time at P0
    time_control_points.push_back(0.5);  // Time at P1
    time_control_points.push_back(1.);  // Time at P2
    cubic_hermite_spline_2Pairs.setTime(time_control_points);
    res1 = cubic_hermite_spline_2Pairs(0.);  // t=0
    ComparePoints(p0, res1, errmsg1, error);
    res1 = cubic_hermite_spline_2Pairs(0.5); // t=0.5
    ComparePoints(p1, res1, errmsg2, error);
    res1 = cubic_hermite_spline_2Pairs(1.);  // t=1
    ComparePoints(p2, res1, errmsg1, error);
    // Test getTime
    try
    {
        cubic_hermite_spline_2Pairs.getTime();
    }
    catch(...)
    {
        error = false;
    }
    if(error)
    {
        std::cout << "Cubic hermite spline test, Error when calling getTime\n";
    }
    // Test derivative : three pairs, time default
    res1 = cubic_hermite_spline_2Pairs.derivate(0.,1);
    ComparePoints(t0, res1, errmsg3, error);
    res1 = cubic_hermite_spline_2Pairs.derivate(0.5,1);
    ComparePoints(t1, res1, errmsg3, error);
    res1 = cubic_hermite_spline_2Pairs.derivate(1.,1);
    ComparePoints(t2, res1, errmsg3, error);
}


void piecewiseCurveTest(bool& error)
{
    // TEST WITH POLYNOMIALS
    std::string errmsg1("in piecewise polynomial curve test, Error While checking value of point on curve : ");
    point_t a(1,1,1); // in [0,1[
    point_t b(2,1,1); // in [1,2[
    point_t c(3,1,1); // in [2,3]
    point_t res;
    t_point_t vec1, vec2, vec3;
    vec1.push_back(a); // x=1, y=1, z=1
    vec2.push_back(b); // x=2, y=1, z=1
    vec3.push_back(c); // x=3, y=1, z=1
    // Create three polynomials of constant value in the interval of definition
    polynomial_t pol1(vec1.begin(), vec1.end(), 0, 1);
    polynomial_t pol2(vec2.begin(), vec2.end(), 1, 2);
    polynomial_t pol3(vec3.begin(), vec3.end(), 2, 3);

    // 1 polynomial in curve
    piecewise_polynomial_curve_t pc(pol1);
    res = pc(0.5);
    ComparePoints(a,res,errmsg1,error);

    // 3 polynomials in curve
    pc.add_curve(pol2);
    pc.add_curve(pol3);

    // Check values on piecewise curve
    // t in [0,1[ -> res=a
    res = pc(0.);
    ComparePoints(a,res,errmsg1,error);
    res = pc(0.5);
    ComparePoints(a,res,errmsg1,error);
    // t in [1,2[ -> res=b
    res = pc(1.0);
    ComparePoints(b,res,errmsg1,error);
    res = pc(1.5);
    ComparePoints(b,res,errmsg1,error);
    // t in [2,3] -> res=c
    res = pc(2.0);
    ComparePoints(c,res,errmsg1,error);
    res = pc(3.0);
    ComparePoints(c,res,errmsg1,error);

    // Create piecewise curve C0 from bezier
    point_t a0(1,2,3);
    point_t b0(2,3,4);
    point_t c0(3,4,5);
    point_t d0(4,5,6);
    std::vector<point_t> params0;
    std::vector<point_t> params1;
    params0.push_back(a0); // bezier between [0,1]
    params0.push_back(b0);
    params0.push_back(c0);
    params0.push_back(d0);
    params1.push_back(d0); // bezier between [1,2]
    params1.push_back(c0); 
    params1.push_back(b0);
    params1.push_back(a0);
    bezier_curve_t bc0(params0.begin(), params0.end(), 0., 1.);
    bezier_curve_t bc1(params1.begin(), params1.end(), 1., 2.);
    piecewise_bezier_curve_t pc_C0(bc0);
    pc_C0.add_curve(bc1);
    // Check value in t=0.5 and t=1.5
    res = pc_C0(0.0);
    ComparePoints(a0,res,errmsg1,error);
    res = pc_C0(1.0);
    ComparePoints(d0,res,errmsg1,error);
    res = pc_C0(2.0);
    ComparePoints(a0,res,errmsg1,error);

    // Create piecewise curve C1 from Hermite
    point_t p0(0.,0.,0.);
    point_t p1(1.,2.,3.);
    point_t p2(4.,4.,4.);
    tangent_t t0(0.5,0.5,0.5);
    tangent_t t1(0.1,0.2,-0.5);
    tangent_t t2(0.1,0.2,0.3);
    std::vector< pair_point_tangent_t > control_points_0;
    control_points_0.push_back(pair_point_tangent_t(p0,t0));
    control_points_0.push_back(pair_point_tangent_t(p1,t1)); // control_points_0 = 1st piece of curve
    std::vector< pair_point_tangent_t > control_points_1;
    control_points_1.push_back(pair_point_tangent_t(p1,t1));
    control_points_1.push_back(pair_point_tangent_t(p2,t2)); // control_points_1 = 2nd piece of curve
    std::vector< double > time_control_points0, time_control_points1;
    time_control_points0.push_back(0.);
    time_control_points0.push_back(1.); // hermite 0 between [0,1]
    time_control_points1.push_back(1.);
    time_control_points1.push_back(3.); // hermite 1 between [1,3]
    cubic_hermite_spline_t chs0(control_points_0.begin(), control_points_0.end(), time_control_points0);
    cubic_hermite_spline_t chs1(control_points_1.begin(), control_points_1.end(), time_control_points1);
    piecewise_cubic_hermite_curve_t pc_C1(chs0);
    pc_C1.add_curve(chs1);


    // Create piecewise curve C2
    point_t a1(0,0,0);
    point_t b1(1,1,1);
    t_point_t veca, vecb;
    // in [0,1[
    veca.push_back(a1);
    veca.push_back(b1); // x=t, y=t, z=t 
    // in [1,2]
    vecb.push_back(b1);
    vecb.push_back(b1); // x=(t-1)+1, y=(t-1)+1, z=(t-1)+1
    polynomial_t pola(veca.begin(), veca.end(), 0, 1);
    polynomial_t polb(vecb.begin(), vecb.end(), 1, 2);
    piecewise_polynomial_curve_t pc_C2(pola);
    pc_C2.add_curve(polb);


    // check C0 continuity
    std::string errmsg2("in piecewise polynomial curve test, Error while checking continuity C0 on ");
    std::string errmsg3("in piecewise polynomial curve test, Error while checking continuity C1 on ");
    std::string errmsg4("in piecewise polynomial curve test, Error while checking continuity C2 on ");
    // not C0
    bool isC0 = pc.is_continuous(0);
    if (isC0)
    {
        std::cout << errmsg2 << " pc " << std::endl;
        error = true;
    }
    // C0
    isC0 = pc_C0.is_continuous(0);
    if (not isC0)
    {
        std::cout << errmsg2 << " pc_C0 " << std::endl;
        error = true;
    }
    // not C1
    bool isC1 = pc_C0.is_continuous(1);
    if (isC1)
    {
        std::cout << errmsg3 << " pc_C0 " << std::endl;
        error = true;
    }
    // C1
    isC1 = pc_C1.is_continuous(1);
    if (not isC1)
    {
        std::cout << errmsg3 << " pc_C1 " << std::endl;
        error = true;
    }
    // not C2
    bool isC2 = pc_C1.is_continuous(2);
    if (isC2)
    {
        std::cout << errmsg4 << " pc_C1 " << std::endl;
        error = true;
    }
    // C2
    isC2 = pc_C2.is_continuous(2);
    if (not isC2)
    {
        std::cout << errmsg4 << " pc_C2 " << std::endl;
        error = true;
    }

    // CONVERT PIECEWISE POLYNOMIAL CURVES TO BEZIER AND HERMITE
    std::string errmsg5("in piecewise polynomial curve test, Error while checking piecewise curve conversion");
    piecewise_bezier_curve_t pc_bezier = pc.convert_piecewise_curve_to_bezier<bezier_curve_t>();
    CompareCurves<piecewise_polynomial_curve_t, piecewise_bezier_curve_t>(pc, pc_bezier, errmsg5, error);
    piecewise_cubic_hermite_curve_t pc_hermite = pc.convert_piecewise_curve_to_cubic_hermite<cubic_hermite_spline_t>();
    CompareCurves<piecewise_polynomial_curve_t, piecewise_cubic_hermite_curve_t>(pc, pc_hermite, errmsg5, error);
    piecewise_polynomial_curve_t pc_polynomial_same = pc.convert_piecewise_curve_to_polynomial<polynomial_t>();
    CompareCurves<piecewise_polynomial_curve_t, piecewise_polynomial_curve_t>(pc, pc_polynomial_same, errmsg5, error);
}

void curveAbcDimDynamicTest(bool& error)
{
    typedef curve_abc<double,double,true> curve_abc_test_t;
    typedef polynomial  <double, double, 3, true> polynomial_test_t;
    typedef exact_cubic <double, double, 3, true> exact_cubic_test_t;
    typedef exact_cubic_test_t::spline_constraints spline_constraints_test_t;
    typedef bezier_curve  <double, double, 3, true> bezier_curve_test_t;
    typedef cubic_hermite_spline <double, double, 3, true> cubic_hermite_spline_test_t;
    curve_abc_test_t * pt_curve_abc;
    // POLYNOMIAL
    point_t a(1,1,1);
    point_t b(2,2,2);
    t_point_t vec;
    vec.push_back(a);
    vec.push_back(b);
    polynomial_test_t pol(vec.begin(), vec.end(), 0, 1);
    try
    {
        pol(0);
        pol(1);
    }
    catch(...)
    {
        error = false;
    }
    // BEZIER
    bezier_curve_test_t bc = bezier_from_curve<bezier_curve_test_t, polynomial_test_t>(pol);
    try
    {
        bc(0);
        bc(1);
    }
    catch(...)
    {
        error = false;
    }
    // CUBIC HERMITE
    cubic_hermite_spline_test_t chs = hermite_from_curve<cubic_hermite_spline_test_t, polynomial_test_t>(pol);
    try
    {
        chs(0);
        chs(1);
    }
    catch(...)
    {
        error = false;
    }
    // EXACT CUBIC : NOT SUPPORTED, problem to fix later
    curves::T_Waypoint waypoints;
    for(double i = 0; i <= 1; i = i + 0.2)
    {
        waypoints.push_back(std::make_pair(i,point_t(i,i,i)));
    }
    std::string errmsg("Error in ExactCubicVelocityConstraintsTest (1); while checking that given wayPoints are crossed (expected / obtained)");
    spline_constraints_test_t constraints;
    constraints.end_vel = point_t(0,0,0);
    constraints.init_vel = point_t(0,0,0);
    constraints.end_acc = point_t(0,0,0);
    constraints.init_acc = point_t(0,0,0);
    exact_cubic_test_t ec(waypoints.begin(), waypoints.end(), constraints);
    try
    {
        ec(0);
        ec(1);
    }
    catch(...)
    {
        error = false;
    }

    // Test with pointer to curve_abc type
    try
    {
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
    }
    catch(...)
    {
        error = false;
    }
    
}

void piecewiseCurveConversionFromDiscretePointsTest(bool& error)
{
    std::string errMsg("piecewiseCurveConversionFromDiscretePointsTest, Error, value on curve is wrong : ");
    point_t p0(0.,0.,0.);
    point_t p1(1.,2.,3.);
    point_t p2(4.,4.,4.);
    point_t p3(10.,10.,10.);
    point_t p_test_0_5 = (p0+p1)/2.0;
    t_point_t points;
    points.push_back(p0);
    points.push_back(p1);
    points.push_back(p2);
    points.push_back(p3);
    double T_min = 1.0;
    double T_max = 3.0;
    double timestep = (T_max-T_min)/double(points.size()-1);
    piecewise_polynomial_curve_t ppc =  piecewise_polynomial_curve_t::
                                        convert_discrete_points_to_polynomial<polynomial_t>(points,T_min,T_max);
    if (!ppc.is_continuous(0))
    {
        std::cout<<"piecewiseCurveConversionFromDiscretePointsTest, Error, piecewise curve is not C0"<<std::endl;
        error = true;
    }
    ComparePoints(p0, ppc(T_min), errMsg, error);
    ComparePoints(p_test_0_5, ppc(T_min+timestep/2.0), errMsg, error);
    ComparePoints(p1, ppc(T_min+timestep), errMsg, error);
    ComparePoints(p2, ppc(T_min+2*timestep), errMsg, error);
    ComparePoints(p3, ppc(T_max), errMsg, error);
}

int main(int /*argc*/, char** /*argv[]*/)
{
    std::cout << "performing tests... \n";
    bool error = false;
    
    CubicFunctionTest(error);
    ExactCubicNoErrorTest(error);
    ExactCubicPointsCrossedTest(error); // checks that given wayPoints are crossed
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
    CubicHermitePairsPositionDerivativeTest(error);
    piecewiseCurveTest(error);
    piecewiseCurveConversionFromDiscretePointsTest(error);
    toPolynomialConversionTest(error);
    cubicConversionTest(error);
    curveAbcDimDynamicTest(error);

    if(error)
    {
        std::cout << "There were some errors\n";
        return -1;
    } else
    {
        std::cout << "no errors found \n";
        return 0;
    }
}
