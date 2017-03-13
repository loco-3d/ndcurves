/**
* \file bezier_curve.h
* \brief class allowing to create a Bezier curve of dimension 1 <= n <= 3.
* \author Steve T.
* \version 0.1
* \date 06/17/2013
*/


#ifndef _CLASS_BEZIERCURVE
#define _CLASS_BEZIERCURVE

#include "curve_abc.h"
#include "bernstein.h"

#include "MathDefs.h"

#include <vector>
#include <stdexcept>

namespace spline
{
/// \class BezierCurve
/// \brief Represents a curve
///
template<typename Time= double, typename Numeric=Time, std::size_t Dim=3, bool Safe=false
, typename Point= Eigen::Matrix<Numeric, Dim, 1> >
struct bezier_curve : public  curve_abc<Time, Numeric, Dim, Safe, Point>
{
	typedef Point 	point_t;
	typedef Time 	time_t;
    typedef Numeric	num_t;
    typedef std::vector<point_t,Eigen::aligned_allocator<point_t> > t_point_t;

/* Constructors - destructors */
	public:
	///\brief Constructor
	///\param PointsBegin, PointsEnd : the points parametering the Bezier curve
	///\TODO : so far size above 3 is ignored
	template<typename In>
	bezier_curve(In PointsBegin, In PointsEnd, const time_t minBound=0, const time_t maxBound=1)
	: minBound_(minBound)
	, maxBound_(maxBound)
	, size_(std::distance(PointsBegin, PointsEnd))
    , bernstein_(spline::makeBernstein<num_t>(size_-1))
	{
        assert(bernstein_.size() == size_);
		In it(PointsBegin);
		if(Safe && (size_<=1 || minBound == maxBound))
		{
            throw std::out_of_range("can't create bezier min bound is higher than max bound"); // TODO
		}
		for(; it != PointsEnd; ++it)
		{
			pts_.push_back(*it);
		}
	}

	///\brief Destructor
	~bezier_curve()
	{
		// NOTHING
	}

	private:
//	bezier_curve(const bezier_curve&);
//  bezier_curve& operator=(const bezier_curve&);
/* Constructors - destructors */

/*Operations*/
	public:
	///  \brief Evaluation of the cubic spline at time t.
	///  \param t : the time when to evaluate the spine
	///  \param return : the value x(t)
    virtual point_t operator()(const time_t t) const
	{
		num_t nT = (t - minBound_) / (maxBound_ - minBound_);
		if(Safe &! (0 <= nT && nT <= 1))
		{
            throw std::out_of_range("can't evaluate bezier curve, out of range"); // TODO
        }
		else
		{
			num_t dt = (1 - nT);
			switch(size_)
			{	
				case 2 :
					return pts_[0] * dt +  nT * pts_[1];
				break;
				case 3 :
					return 	pts_[0] * dt * dt 
                       				+ 2 * pts_[1] * nT * dt
						+ pts_[2] * nT * nT;
				break;
                case 4 :
					return 	pts_[0] * dt * dt * dt
						+ 3 * pts_[1] * nT * dt * dt 
						+ 3 * pts_[2] * nT * nT * dt 
						+ pts_[3] * nT * nT *nT;
                default :
                    return evalBernstein(dt);
				break;
			}
		}
	}

    ///
    /// \brief Evaluates all Bernstein polynomes for a certain degree
    ///
    point_t evalBernstein(const Numeric u) const
    {
        point_t res = point_t::Zero();
        typename t_point_t::const_iterator pts_it = pts_.begin();
        for(typename std::vector<Bern<Numeric> >::const_iterator cit = bernstein_.begin();
            cit !=bernstein_.end(); ++cit, ++pts_it)
        {
            res += cit->operator()(u) * (*pts_it);
        }
        return res;
    }

/*Operations*/

/*Helpers*/
	virtual time_t min() const{return minBound_;}
	virtual time_t max() const{return maxBound_;}
/*Helpers*/

	public:
    const time_t minBound_, maxBound_;
    const std::size_t size_;
    const std::vector<Bern<Numeric> > bernstein_;
	
    private:
    t_point_t  pts_;

    //storing bernstein polynoms, even in low dimension
};
}
#endif //_CLASS_BEZIERCURVE

