#include "API/CubicFunction.h"
#include "API/SplineVisitor.h"

using namespace spline;

CubicFunction::CubicFunction(const Vector3& a, const Vector3& b, const Vector3& c, const Vector3& d, const Real tBegin, const Real tEnd)
	: a_(a), b_(b), c_(c), d_(d), tBegin_(tBegin), tEnd_(tEnd)
{
	// NOTHING
}

CubicFunction::~CubicFunction()
{
	// NOTHING
}

bool CubicFunction::Evaluate(const Real t, Vector3& result) const
{
	if(tBegin_ <= t && t <= tEnd_)
	{
		Real dt = (t - tBegin_);
		result =  a_ + b_ * dt +  c_ * dt * dt + d_ * dt * dt * dt;
		return true;
	}
	else // t out of bounds
	{
		return false;
	}
}

void CubicFunction::Accept(SplineVisitor& visitor, Real dt) const
{
	for(Real ti = tBegin_; ti <= tEnd_; ti = ti + dt)
	{
		Vector3 res; Evaluate(ti, res);
		visitor.Visit(ti, res);
	}
}
