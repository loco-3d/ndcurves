
#include "CubicFunction.h"
#include "ExactCubic.h"
#include "SplineVisitor.h"

#include <string>
#include <iostream>
#include <cmath>

using namespace std;

namespace spline
{
bool QuasiEqual(const Real a, const Real b, const float margin)
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

bool operator ==(const Vector3& a, const Vector3& b)
{
	return QuasiEqual(a.x(), b.x(), margin) && QuasiEqual(a.y(), b.y(), margin) && QuasiEqual(a.z(), b.z(), margin);
}
} // namespace spline

using namespace spline;

ostream& operator<<(ostream& os, const Vector3& pt)
{
    os << "(" << pt.x() << ", " << pt.y() << ", " << pt.z() << ")";
    return os;
}

void ComparePoints(const Vector3& pt1, const Vector3& pt2, const std::string& errmsg, bool& error)
{
	if(!(pt1 == pt2))
	{
		error = true;
		std::cout << errmsg << pt1 << " ; " << pt2 << "\n";
	}
}

/*Cubic Function tests*/

void CubicFunctionTest(bool& error)
{
	std::string errMsg("In test CubicFunctionTest ; unexpected result for x ");
	Vector3 a(1,2,3);
	Vector3 b(2,3,4);
	Vector3 c(3,4,5);
	Vector3 d(3,6,7);

	CubicFunction cf(a, b, c, d, 0, 1);
	Vector3 res1;
	cf.Evaluate(0, res1); 
	Vector3 x0(1,2,3);
	ComparePoints(x0, res1, errMsg + "(0) ", error);
	
	Vector3 x1(9,15,19);
	cf.Evaluate(1, res1); 
	ComparePoints(x1, res1, errMsg + "(1) ", error);
	
	Vector3 x2(3.125,5.25,7.125);
	cf.Evaluate(0.5, res1); 
	ComparePoints(x2, res1, errMsg + "(0.5) ", error);

	CubicFunction cf2(a, b, c, d, 0.5, 1);
	cf2.Evaluate(0.5, res1); 
	ComparePoints(x0, res1, errMsg + "x3 ", error);
	if(cf2.Evaluate(0.4, res1))
	{
		error = true;
		std::cout << "Evaluation of cubic cf2 error, 0.4 should be an out of range value\n";
	}
	if(cf2.Evaluate(1.1, res1))
	{
		error = true;
		std::cout << "Evaluation of cubic cf2 error, 1.1 should be an out of range value\n";
	}
}

/*Exact Cubic Function tests*/
void ExactCubicNoErrorTest(bool& error)
{
	spline::T_Waypoint waypoints;
	for(Real i = 0; i <= 1; i = i + 0.2)
	{
		waypoints.push_back(std::make_pair(i,Vector3(i,i,i)));
	}
	ExactCubic exactCubic(waypoints);
	Vector3 res1;
	if(!exactCubic.Evaluate(0, res1))
	{
		error = true;
		std::cout << "Evaluation of exactCubic error, 0 should be in range value\n";
	}
	if(!exactCubic.Evaluate(1, res1))
	{
		error = true;
		std::cout << "Evaluation of exactCubic error, 1 should be in range value\n";
	}
	if(exactCubic.Evaluate(1.2, res1))
	{
		error = true;
		std::cout << "Evaluation of exactCubic error, 1.2 should be an out of range value\n";
	}
}

void ExactCubicPointsCrossedTest(bool& error)
{
	spline::T_Waypoint waypoints;
	for(Real i = 0; i <= 1; i = i + 0.2)
	{
		waypoints.push_back(std::make_pair(i,Vector3(i,i,i)));
	}
	ExactCubic exactCubic(waypoints);
	Vector3 res1;
	for(Real i = 0; i <= 1; i = i + 0.2)
	{
		exactCubic.Evaluate(i, res1);
		std::string errmsg("Error While checking that given wayPoints are crossed (expected / obtained)");
		ComparePoints(Vector3(i,i,i), res1, errmsg, error);
	}
}

/*Cubic Visitor tests*/
#include <vector>

#include<Eigen/StdVector>

namespace spline
{
	typedef std::vector<Vector3,Eigen::aligned_allocator<Vector3> > T_Vector;
	typedef T_Vector::const_iterator CIT_Vector;
	struct SplineVisitorTest : public SplineVisitor
	{
		SplineVisitorTest()
			: previousTime_(-1)
		{
			// NOTHING
		}

		~SplineVisitorTest()
		{
			// NOTHING
		}

		virtual void Visit(const Real time, const Vector3& value)
		{
			if(previousTime_ >= time)
			{
				std::cout << "Error : Visit method called with non sequential time values" << std::endl;
			}
			previousTime_ = time;
			values_.push_back(value);
		}

		Real previousTime_;
		T_Vector values_;
	};
}

void SplineVisitorTestFunction(bool& error)
{
	spline::T_Waypoint waypoints;
	for(Real i = 0; i <= 1; i = i + 0.2)
	{
		waypoints.push_back(std::make_pair(i,Vector3(i,i,i)));
	}
	ExactCubic exactCubic(waypoints);
	Vector3 res1;
	SplineVisitorTest visitor;
	exactCubic.Accept(visitor, 0.2);
	CIT_Vector it = visitor.values_.begin();
	for(Real i = 0; i <= 1; i = i + 0.2)
	{
		assert(it != visitor.values_.end());
		exactCubic.Evaluate(i, res1);
		std::string errmsg("Error While testing SplineVisitor at timestep (expected / obtained)");
		ComparePoints(*it, res1, errmsg, error);
		++it;
	}
	error = error & exactCubic.Evaluate(0.9, res1);
	std::string errmsg("Error While testing SplineVisitor at timestep (expected / obtained)");
	ComparePoints(Vector3(0.923, 0.923 ,0.923), res1, errmsg, error);
}

int main(int argc, char *argv[])
{
	std::cout << "performing tests... \n";
	bool error = false;
	CubicFunctionTest(error);
	ExactCubicNoErrorTest(error);
	ExactCubicPointsCrossedTest(error); // checks that given wayPoints are crossed
	SplineVisitorTestFunction(error); // checks that given wayPoints are crossed
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
