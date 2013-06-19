#include "API/ExactCubic.h"
#include "API/CubicFunction.h"
#include "API/SplineVisitor.h"

#include <vector>

namespace spline
{
typedef std::pair<Real, CubicFunction*>	SubSpline;
typedef std::vector<SubSpline>			T_SubSpline;
typedef T_SubSpline::iterator			IT_SubSpline;
typedef T_SubSpline::const_iterator		CIT_SubSpline;

struct CubicPImpl{
	CubicPImpl(Real maxTime)
		: maxTime_(maxTime)
	{
		// NOTHING
	}

	~CubicPImpl()
	{
		for(IT_SubSpline it = subSplines_.begin(); it != subSplines_.end(); ++ it)
		{
			delete it->second;
		}
	}

	T_SubSpline subSplines_;
	Real maxTime_; // max bound on which spline is defined.
};
} // namespace spline

using namespace spline;

ExactCubic::ExactCubic(const T_Waypoint& waypoints)
{
	size_t size = waypoints.size();
	assert(size > 2);
	pImpl_.reset(new CubicPImpl(waypoints[size-1].first)); // setting max boundary
	// refer to the paper to understand all this.
	MatrixX h1 = MatrixX::Zero(size, size);
	MatrixX h2 = MatrixX::Zero(size, size);
	MatrixX h3 = MatrixX::Zero(size, size);
	MatrixX h4 = MatrixX::Zero(size, size);
	MatrixX h5 = MatrixX::Zero(size, size);
	MatrixX h6 = MatrixX::Zero(size, size);

	MatrixX a =  MatrixX::Zero(size, 3);
	MatrixX b =  MatrixX::Zero(size, 3);
	MatrixX c =  MatrixX::Zero(size, 3);
	MatrixX d =  MatrixX::Zero(size, 3);
	MatrixX x =  MatrixX::Zero(size, 3);

	Real dTi, dTi_1, dTi_sqr, dTi_1sqr, dTi_cube;
	Real t_previous = waypoints[0].first;

	// I think using numbers instead of iterators will be clearer on this case
	for(int i=0; i< size - 1; ++i)
	{
		dTi      = waypoints[i + 1].first - waypoints[i].first;
		dTi_sqr  = dTi * dTi;
		dTi_cube = dTi_sqr * dTi;
		
		assert(dTi > 0); //make sure of the ascendant order

		// filling matrices values
		h3(i,i)   = -3 / dTi_sqr;
		h3(i,i+1) =  3 / dTi_sqr;
				
		h4(i,i)   = -2 / dTi;
		h4(i,i+1) = -1 / dTi;

		h5(i,i)   =  2 / dTi_cube;
		h5(i,i+1) = -2 / dTi_cube;
		
		h6(i,i)   = 1 / dTi_sqr;
		h6(i,i+1) = 1 / dTi_sqr;

		// we stop one step earlier for matrices h1 and h2
		if(i + 2 < size)
		{
			// this can be optimized but let's focus on clarity as long as not needed
			dTi_1	  = waypoints[i + 2].first - waypoints[i + 1].first;
			dTi_1sqr  = dTi_1 * dTi_1;

			h1(i+1, i)   = 2 / dTi;
			h1(i+1, i+1) = 4 / dTi + 4 / dTi_1;
			h1(i+1, i+2) = 2 / dTi_1;

			h2(i+1, i)   = - 6 / dTi_sqr;
			h2(i+1, i+1) =   (6 / dTi_1sqr) - (6 / dTi_sqr);
			h2(i+1, i+2) =   6 / dTi_1sqr;
		}

		//computing x = [x0* x1* ... xT*]^T
		x.row(i) =  waypoints[i].second.transpose();
	}
	// adding last x
	x.row(size-1) =  waypoints[size-1] .second.transpose();

	// now we just have to apply the linear relations:

	a = x;
	// should I pseudo inverse this?
	PseudoInverse(h1);
	b = h1 * h2 * x; //h1 * b = h2 * x => b = (h1)^-1 * h2 * x
	c = h3 * x + h4 * b;
	d = h5 * x + h6 * b;

	// Now each line i of the a, b, c and d matrices correponds to the coefficient of spline i
	// Let's create those, and initialiaze our pImpl_
	for(int i=0; i< size - 1; ++i) // last spline never evaluated, that's ok because xi(ti_1) = x_i_1* 
	{
		CubicFunction* subSpline = new CubicFunction(a.row(i), b.row(i), c.row(i), d.row(i), waypoints[i].first, waypoints[i+1].first);
		pImpl_->subSplines_.push_back(std::make_pair(waypoints[i].first, subSpline));
	}
	
}

ExactCubic::~ExactCubic()
{
	// NOTHING
}


bool ExactCubic::Evaluate(const Real t, Vector3& value) const
{	
	CIT_SubSpline it = pImpl_->subSplines_.begin();
	CIT_SubSpline it2 = it;
	++it2;
	for(; it2!= pImpl_->subSplines_.end(); ++it2)
	{
		// t out of min bound
		if(t < it->first)
		{
			return false;
		}
		else if(t <= it2->first)
		{
			return it->second->Evaluate(t, value);
		}
		++it;
	}
	// handling max Bound
	if(t <= pImpl_->maxTime_)
	{
		return it->second->Evaluate(t, value);
	}
	else
	{
		return false; // t out of max bound
	}
}

void ExactCubic::Accept(SplineVisitor& visitor, Real dt) const
{
	for(Real ti = pImpl_->subSplines_[0].first; ti <= pImpl_->maxTime_; ti = ti + dt)
	{
		Vector3 res; Evaluate(ti, res);
		visitor.Visit(ti, res);
	}
}
