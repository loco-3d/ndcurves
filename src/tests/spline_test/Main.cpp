
#include "spline/exact_cubic.h"
#include "spline/bezier_curve.h"
#include "spline/spline_curve.h"

#include <string>
#include <iostream>
#include <cmath>

using namespace std;

namespace spline
{
typedef Eigen::Vector3d point_t;
typedef std::vector<point_t,Eigen::aligned_allocator<point_t> >  t_point_t;
typedef spline_curve  <double, double, 3, 3, true, point_t, t_point_t> cubic_function_t;
typedef exact_cubic   <double, double, 3, true, point_t> exact_cubic_t;
typedef bezier_curve  <double, double, 3, true, point_t> bezier_curve_t;
typedef std::pair<double, point_t> Waypoint;
typedef std::vector<Waypoint> T_Waypoint;


typedef Eigen::Matrix<double,1,1> point_one;
typedef spline_curve<double, double, 1, 3, true, point_one> cubic_function_one;
typedef exact_cubic   <double, double, 1, true, point_one> exact_cubic_one;
typedef std::pair<double, point_one> WaypointOne;
typedef std::vector<WaypointOne> T_WaypointOne;

bool QuasiEqual(const double a, const double b, const float margin)
{
	if ((a <= 0 && b <= 0) || (a >= 0 && b>= 0))
	{
		return (abs(a-b)) <= margin;
	}
	else
	{
		return abs(a) + abs(b) <= margin;
	}
}

const float margin = 0.01f;

bool operator ==(const point_t& a, const point_t& b)
{
	return QuasiEqual(a.x(), b.x(), margin) && QuasiEqual(a.y(), b.y(), margin) && QuasiEqual(a.z(), b.z(), margin);
}

} // namespace spline

using namespace spline;

ostream& operator<<(ostream& os, const point_t& pt)
{
    os << "(" << pt.x() << ", " << pt.y() << ", " << pt.z() << ")";
    return os;
}

void ComparePoints(const point_t& pt1, const point_t& pt2, const std::string& errmsg, bool& error)
{
	if(!(pt1 == pt2))
	{
		error = true;
        std::cout << errmsg << pt1 << " ; " << pt2 << std::endl;
	}
}

void ComparePoints(const point_one& pt1, const point_one& pt2, const std::string& errmsg, bool& error)
{
	if(!(pt1 == pt2))
	{
		error = true;
        std::cout << errmsg << pt1 << " ; " << pt2 <<  std::endl;
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
    cubic_function_t cf(vec.begin(), vec.end(), 0, 1);
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
    cubic_function_t cf2(vec, 0.5, 1);
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
	if(cf.max() != 1)
	{
		error = true;
		std::cout << "Evaluation of exactCubic error, MaxBound should be equal to 1\n";
	}
	if(cf.min() != 0)
	{
		error = true;
		std::cout << "Evaluation of exactCubic error, MinBound should be equal to 1\n";
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
	params.push_back(b);

	// 2d curve
	bezier_curve_t cf(params.begin(), params.end());
	point_t res1;
	res1 = cf(0); 
	point_t x20 = a ;
    ComparePoints(x20, res1, errMsg + "2(0) ", error);
	
	point_t x21 = b;
	res1 = cf(1); 
    ComparePoints(x21, res1, errMsg + "2(1) ", error);
	
	//3d curve
	params.push_back(c);
	bezier_curve_t cf3(params.begin(), params.end());
	res1 = cf3(0); 
    ComparePoints(a, res1, errMsg + "3(0) ", error);

	res1 = cf3(1); 
    ComparePoints(c, res1, errMsg + "3(1) ", error);

	//4d curve
	params.push_back(d);
	bezier_curve_t cf4(params.begin(), params.end(), 0.4, 2);
	res1 = cf4(0.4); 
    ComparePoints(a, res1, errMsg + "3(0) ", error);
	
	res1 = cf4(2); 
    ComparePoints(d, res1, errMsg + "3(1) ", error);

	try
	{
		cf(-0.4);
	}
	catch(...)
	{
		error = false;
	}
	if(error)
	{
		std::cout << "Evaluation of bezier cf error, -0.4 should be an out of range value\n";
	}
	error = true;	
	try
	{
		cf(1.1);
	}
	catch(...)
	{
		error = false;
	}
	if(error)
	{
		std::cout << "Evaluation of bezier cf error, 1.1 should be an out of range value\n";
	}
	if(cf.max() != 1)
	{
		error = true;
		std::cout << "Evaluation of exactCubic error, MaxBound should be equal to 1\n";
	}
	if(cf.min() != 0)
	{
		error = true;
		std::cout << "Evaluation of exactCubic error, MinBound should be equal to 1\n";
	}
}

/*Exact Cubic Function tests*/
void ExactCubicNoErrorTest(bool& error)
{
	spline::T_Waypoint waypoints;
	for(double i = 0; i <= 1; i = i + 0.2)
	{
		waypoints.push_back(std::make_pair(i,point_t(i,i,i)));
	}
	exact_cubic_t exactCubic(waypoints.begin(), waypoints.end());
	point_t res1;
	try
	{
		exactCubic(0);
		exactCubic(1);
	}
	catch(...)
	{
		error = true;
		std::cout << "Evaluation of ExactCubicNoErrorTest error\n";
	}
	error = true;
	try
	{
		exactCubic(1.2);
	}
	catch(...)
	{
		error = false;
	}
	if(error)
	{
		std::cout << "Evaluation of exactCubic cf error, 1.2 should be an out of range value\n";
	}
	if(exactCubic.max() != 1)
	{
		error = true;
		std::cout << "Evaluation of exactCubic error, MaxBound should be equal to 1\n";
	}
	if(exactCubic.min() != 0)
	{
		error = true;
		std::cout << "Evaluation of exactCubic error, MinBound should be equal to 1\n";
	}
}


/*Exact Cubic Function tests*/
void ExactCubicTwoPointsTest(bool& error)
{
	spline::T_Waypoint waypoints;
	for(double i = 0; i < 2; ++i)
	{
		waypoints.push_back(std::make_pair(i,point_t(i,i,i)));
	}
	exact_cubic_t exactCubic(waypoints.begin(), waypoints.end());
	
	point_t res1 = exactCubic(0);
	std::string errmsg("in ExactCubic 2 points Error While checking that given wayPoints  are crossed (expected / obtained)");
	ComparePoints(point_t(0,0,0), res1, errmsg, error);
		
	res1 = exactCubic(1);
	ComparePoints(point_t(1,1,1), res1, errmsg, error);
}

void ExactCubicOneDimTest(bool& error)
{
	spline::T_WaypointOne waypoints;
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

void ExactCubicPointsCrossedTest(bool& error)
{
	spline::T_Waypoint waypoints;
	for(double i = 0; i <= 1; i = i + 0.2)
	{
		waypoints.push_back(std::make_pair(i,point_t(i,i,i)));
	}
	exact_cubic_t exactCubic(waypoints.begin(), waypoints.end());
	point_t res1;
	for(double i = 0; i <= 1; i = i + 0.2)
	{
		res1 = exactCubic(i);
		std::string errmsg("Error While checking that given wayPoints are crossed (expected / obtained)");
		ComparePoints(point_t(i,i,i), res1, errmsg, error);
	}
}

int main(int argc, char *argv[])
{
	std::cout << "performing tests... \n";
	bool error = false;
	CubicFunctionTest(error);
	ExactCubicNoErrorTest(error);
	ExactCubicPointsCrossedTest(error); // checks that given wayPoints are crossed
	ExactCubicTwoPointsTest(error);
	ExactCubicOneDimTest(error);
	//BezierCurveTest(error);
	if(error)
	{
		std::cout << "There were some errors\n";
		return -1;
	}
	else
	{
		std::cout << "no errors found \n";
		return 0;
	}
}

