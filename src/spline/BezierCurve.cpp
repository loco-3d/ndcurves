#include "API/BezierCurve.h"

using namespace spline;

BezierCurve::BezierCurve(const T_Vector&  points, const Real minBound, const Real maxBound)
	: size_  ((int)points.size())
	, pts_(points)
	, minBound_(minBound)
	, maxBound_(maxBound)
{
	assert(size_>1 && minBound < maxBound);
	// NOTHING
}

BezierCurve::~BezierCurve()
{
	// NOTHING
}

bool BezierCurve::Evaluate(const Real t, Vector3& result) const
{
	Real nT = (t - minBound_) / (maxBound_ - minBound_);
	if(0 <= nT && nT <= 1)
	{
		Real dt = (1 - nT);
		switch(size_)
		{	
			case 2 :
				result =  pts_[0] * dt +  nT * pts_[1];
			break;
			case 3 :
				result =  pts_[0] * dt * dt +  2 * pts_[1] * nT * dt  + pts_[2] * nT * nT ;
			break;
			default :
				result =  pts_[0] * dt * dt * dt +  3 * pts_[1] * nT * dt * dt + 3 * pts_[2] * nT * nT * dt + pts_[3] * nT * nT *nT ;
			break;
		}
		return true;
	}
	else // t out of bounds
	{
		return false;
	}
}

Real BezierCurve::MinBound() const
{
	return minBound_;
}

Real BezierCurve::MaxBound() const
{
	return maxBound_;
}