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
#include "curve_constraint.h"

#include "MathDefs.h"

#include <vector>
#include <stdexcept>

#include <iostream>

namespace spline
{
/// \class BezierCurve
/// \brief Represents a Bezier curve of arbitrary dimension and order.
/// For degree lesser than 4, the evaluation is analitycal.Otherwise
/// the bernstein polynoms are used to evaluate the spline at a given location.
///
template<typename Time= double, typename Numeric=Time, std::size_t Dim=3, bool Safe=false
, typename Point= Eigen::Matrix<Numeric, Dim, 1> >
struct bezier_curve : public curve_abc<Time, Numeric, Dim, Safe, Point>
{
	typedef Point 	point_t;
	typedef Time 	time_t;
    typedef Numeric	num_t;
    typedef curve_constraints<point_t> curve_constraints_t;
    typedef std::vector<point_t,Eigen::aligned_allocator<point_t> > t_point_t;
    typedef typename t_point_t::const_iterator cit_point_t;
    typedef bezier_curve<Time, Numeric, Dim, Safe, Point > bezier_curve_t;

/* Constructors - destructors */
	public:
	///\brief Constructor
	///\param PointsBegin, PointsEnd : the points parametering the Bezier curve
    ///
	template<typename In>
	bezier_curve(In PointsBegin, In PointsEnd, const time_t minBound=0, const time_t maxBound=1)
	: minBound_(minBound)
	, maxBound_(maxBound)
	, size_(std::distance(PointsBegin, PointsEnd))
    , degree_(size_-1)
    , bernstein_(spline::makeBernstein<num_t>(degree_))
    {
        assert(bernstein_.size() == size_);
		In it(PointsBegin);
        if(Safe && (size_<1 || minBound >= maxBound))
            throw std::out_of_range("can't create bezier min bound is higher than max bound"); // TODO
        for(; it != PointsEnd; ++it)
            pts_.push_back(*it);
    }


    ///\brief Constructor
    /// This constructor will add 4 points (2 after the first one, 2 before the last one)
    /// to ensure that velocity and acceleration constraints are respected
    ///\param PointsBegin, PointsEnd : the points parametering the Bezier curve
    ///\param constraints : constraints applying on start / end velocities and acceleration
    ///
    template<typename In>
    bezier_curve(In PointsBegin, In PointsEnd, const curve_constraints_t& constraints, const time_t minBound=0, const time_t maxBound=1)
    : minBound_(minBound)
    , maxBound_(maxBound)
    , size_(std::distance(PointsBegin, PointsEnd)+4)
    , degree_(size_-1)
    , bernstein_(spline::makeBernstein<num_t>(degree_))
    {
        if(Safe && (size_<1 || minBound >= maxBound))
            throw std::out_of_range("can't create bezier min bound is higher than max bound");
        t_point_t updatedList = add_constraints<In>(PointsBegin, PointsEnd, constraints);
        for(cit_point_t cit = updatedList.begin(); cit != updatedList.end(); ++cit)
            pts_.push_back(*cit);
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
                case 1 :
                    return pts_[0];
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
                    return evalHorner(nT);
				break;
			}
		}
	}

    ///  \brief Computes the derivative curve at order N.
    ///  \param order : order of the derivative
    ///  \param return : the value x(t)
    bezier_curve_t compute_derivate(const std::size_t order) const
    {
        if(order == 0) return *this;
        t_point_t derived_wp;
        for(typename t_point_t::const_iterator pit =  pts_.begin(); pit != pts_.end()-1; ++pit)
            derived_wp.push_back(degree_ * (*(pit+1) - (*pit)));
        if(derived_wp.empty())
            derived_wp.push_back(point_t::Zero());
        bezier_curve_t deriv(derived_wp.begin(), derived_wp.end(),minBound_,maxBound_);
        return deriv.compute_derivate(order-1);
    }

    ///  \brief Computes the primitive of the curve at order N.
    ///  \param constant : value of the primitive at t = 0
    ///  \param return : the curve x_1(t) such that d/dt(x_1(t)) = x_1(t)
    bezier_curve_t compute_primitive(const std::size_t order) const
    {
        if(order == 0) return *this;
        num_t new_degree = (num_t)(degree_+1);
        t_point_t n_wp;
        point_t current_sum =  point_t::Zero();
        // recomputing waypoints q_i from derivative waypoints p_i. q_0 is the given constant.
        // then q_i = (sum( j = 0 -> j = i-1) p_j) /n+1
        n_wp.push_back(current_sum);
        for(typename t_point_t::const_iterator pit =  pts_.begin(); pit != pts_.end(); ++pit)
        {
            current_sum += *pit;
            n_wp.push_back(current_sum / new_degree);
        }
        bezier_curve_t integ(n_wp.begin(), n_wp.end(),minBound_,maxBound_);
        return integ.compute_primitive(order-1);
    }

    ///  \brief Evaluates the derivative at order N of the curve.
    ///  If the derivative is to be evaluated several times, it is
    ///  rather recommended to compute the derivative curve using compute_derivate
    ///  \param order : order of the derivative
    ///  \param t : the time when to evaluate the spine
    ///  \param return : the value x(t)
    virtual point_t derivate(const time_t t, const std::size_t order) const
    {
        bezier_curve_t deriv =compute_derivate(order);
        return deriv(t);
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
            res += cit->operator()(u) * (*pts_it);
        return res;
    }


    ///
    /// \brief Evaluates all Bernstein polynomes for a certain degree using horner's scheme
    ///
    point_t evalHorner(const Numeric t) const
    {
        typename t_point_t::const_iterator pts_it = pts_.begin();
        Numeric u, bc, tn;
        u = 1.0 - t;
        bc = 1;
        tn = 1;
        point_t tmp =(*pts_it)*u; ++pts_it;
        for(int i=1; i<degree_; i++, ++pts_it)
        {
            tn = tn*t;
            bc = bc*(degree_-i+1)/i;
            tmp = (tmp + tn*bc*(*pts_it))*u;
        }
        return (tmp + tn*t*(*pts_it));
    }

    const t_point_t& waypoints() const {return pts_;}

    private:
    template<typename In>
    t_point_t add_constraints(In PointsBegin, In PointsEnd, const curve_constraints_t& constraints)
    {
        t_point_t res;
        point_t P0, P1, P2, P_n_2, P_n_1, PN;
        P0 = *PointsBegin; PN = *(PointsEnd-1);
        P1    = P0+ constraints.init_vel / degree_;
        P_n_1 = PN -constraints.end_vel  / degree_;
        P2    = constraints.init_acc / (degree_ * (degree_-1)) + 2* P1    - P0;
        P_n_2 = constraints.end_acc  / (degree_ * (degree_-1)) + 2* P_n_1 - PN;

        res.push_back(P0);
        res.push_back(P1);
        res.push_back(P2);

        for(In it = PointsBegin+1; it != PointsEnd-1; ++it)
            res.push_back(*it);

        res.push_back(P_n_2);
        res.push_back(P_n_1);
        res.push_back(PN);
        return res;
    }

/*Operations*/

/*Helpers*/
    public:
	virtual time_t min() const{return minBound_;}
	virtual time_t max() const{return maxBound_;}
/*Helpers*/

	public:
    const time_t minBound_, maxBound_;
    const std::size_t size_;
    const std::size_t degree_;
    const std::vector<Bern<Numeric> > bernstein_;
	
    private:
    t_point_t  pts_;

    //storing bernstein polynoms, even in low dimension
};
}
#endif //_CLASS_BEZIERCURVE

