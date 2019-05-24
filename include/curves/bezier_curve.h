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

namespace curves
{
/// \class BezierCurve.
/// \brief Represents a Bezier curve of arbitrary dimension and order.
/// For degree lesser than 4, the evaluation is analitycal. Otherwise
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

    /// \brief Constructor.
    /// Given the first and last point of a control points set, create the bezier curve.
    /// \param PointsBegin  : an iterator pointing to the first element of a control points container.
    /// \param PointsEnd    : an iterator pointing to the last element of a control points container.
    ///
	template<typename In>
    bezier_curve(In PointsBegin, In PointsEnd)
    : T_min_(0.)
    , T_max_(1.)
    , mult_T_(1.)
	, size_(std::distance(PointsBegin, PointsEnd))
    , degree_(size_-1)
    , bernstein_(curves::makeBernstein<num_t>((unsigned int)degree_))
    {
        assert(bernstein_.size() == size_);
		In it(PointsBegin);
        if(Safe && (size_<1 || T_max_ <= T_min_)) 
        {
            throw std::out_of_range("can't create bezier min bound is higher than max bound");
        }
        for(; it != PointsEnd; ++it)
        {
            pts_.push_back(*it);
        }
    }

    /// \brief Constructor.
    /// Given the first and last point of a control points set, create the bezier curve.
    /// \param PointsBegin   : an iterator pointing to the first element of a control point container.
    /// \param PointsEnd     : an iterator pointing to the last element of a control point container.
    /// \param T             : upper bound of curve parameter which is between \f$[0;T]\f$ (default \f$[0;1]\f$).
    ///
    template<typename In>
    bezier_curve(In PointsBegin, In PointsEnd, const time_t T_min, const time_t T_max)
    : T_min_(T_min)
    , T_max_(T_max)
    , mult_T_(1.)
    , size_(std::distance(PointsBegin, PointsEnd))
    , degree_(size_-1)
    , bernstein_(curves::makeBernstein<num_t>((unsigned int)degree_))
    {
        assert(bernstein_.size() == size_);
        In it(PointsBegin);
        if(Safe && (size_<1 || T_max_ <= T_min_))
        {
            throw std::out_of_range("can't create bezier min bound is higher than max bound"); // TODO
        }
        for(; it != PointsEnd; ++it)
        {
            pts_.push_back(*it);
        }
    }



    /// \brief Constructor.
    /// Given the first and last point of a control points set, create the bezier curve.
    /// \param PointsBegin   : an iterator pointing to the first element of a control point container.
    /// \param PointsEnd     : an iterator pointing to the last element of a control point container.
    /// \param T             : upper bound of time which is between \f$[0;T]\f$ (default \f$[0;1]\f$).
    /// \param mult_T        : ... (default value is 1.0).
    ///
    template<typename In>
    bezier_curve(In PointsBegin, In PointsEnd, const time_t T_min, const time_t T_max, const time_t mult_T)
    : T_min_(T_min)
    , T_max_(T_max)
    , mult_T_(mult_T)
    , size_(std::distance(PointsBegin, PointsEnd))
    , degree_(size_-1)
    , bernstein_(curves::makeBernstein<num_t>((unsigned int)degree_))
    {
        assert(bernstein_.size() == size_);
        In it(PointsBegin);
        if(Safe && (size_<1 || T_max_ <= T_min_))
        {
            throw std::out_of_range("can't create bezier min bound is higher than max bound"); // TODO
        }
        for(; it != PointsEnd; ++it)
        {
            pts_.push_back(*it);
        }
    }

    /// \brief Constructor
    /// This constructor will add 4 points (2 after the first one, 2 before the last one)
    /// to ensure that velocity and acceleration constraints are respected.
    /// \param PointsBegin   : an iterator pointing to the first element of a control point container.
    /// \param PointsEnd     : an iterator pointing to the last element of a control point container.
    /// \param constraints : constraints applying on start / end velocities and acceleration.
    ///
    template<typename In>
    bezier_curve(In PointsBegin, In PointsEnd, const curve_constraints_t& constraints, 
                const time_t T_min=0., const time_t T_max=1.)
    : T_min_(T_min)
    , T_max_(T_max)
    , mult_T_(1.)
    , size_(std::distance(PointsBegin, PointsEnd)+4)
    , degree_(size_-1)
    , bernstein_(curves::makeBernstein<num_t>((unsigned int)degree_))
    {
        if(Safe && (size_<1 || T_max_ <= T_min_))
        {
            throw std::out_of_range("can't create bezier min bound is higher than max bound");
        }
        t_point_t updatedList = add_constraints<In>(PointsBegin, PointsEnd, constraints);
        for(cit_point_t cit = updatedList.begin(); cit != updatedList.end(); ++cit)
        {
            pts_.push_back(*cit);
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
	///  \brief Evaluation of the bezier curve at time t.
	///  \param t : time when to evaluate the curve.
	///  \return \f$x(t)\f$ point corresponding on curve at time t.
    virtual point_t operator()(const time_t t) const
    {
        if(Safe &! (T_min_ <= t && t <= T_max_))
        {
            throw std::out_of_range("can't evaluate bezier curve, out of range"); // TODO
        }
        if (size_ == 1)
        {
          return mult_T_*pts_[0];
        }else
        {
          return evalHorner(t);
        }
	}

    ///  \brief Compute the derived curve at order N.
    ///  Computes the derivative order N, \f$\frac{d^Nx(t)}{dt^N}\f$ of bezier curve of parametric equation x(t).
    ///  \param order : order of derivative.
    ///  \return \f$\frac{d^Nx(t)}{dt^N}\f$ derivative order N of the curve.
    bezier_curve_t compute_derivate(const std::size_t order) const
    {
        if(order == 0) 
        {
            return *this;
        }
        t_point_t derived_wp;
        for(typename t_point_t::const_iterator pit =  pts_.begin(); pit != pts_.end()-1; ++pit)
            derived_wp.push_back((num_t)degree_ * (*(pit+1) - (*pit)));
        if(derived_wp.empty())
            derived_wp.push_back(point_t::Zero(Dim));
        bezier_curve_t deriv(derived_wp.begin(), derived_wp.end(),T_min_, T_max_, mult_T_ * (1./(T_max_-T_min_)) );
        return deriv.compute_derivate(order-1);
    }

    ///  \brief Compute the primitive of the curve at order N.
    ///  Computes the primitive at order N of bezier curve of parametric equation \f$x(t)\f$. <br>
    ///  At order \f$N=1\f$, the primitve \f$X(t)\f$ of \f$x(t)\f$ is such as \f$\frac{dX(t)}{dt} = x(t)\f$.
    ///  \param order : order of the primitive.
    ///  \return primitive at order N of x(t).
    bezier_curve_t compute_primitive(const std::size_t order) const
    {
        if(order == 0) 
        {
            return *this;
        }
        num_t new_degree = (num_t)(degree_+1);
        t_point_t n_wp;
        point_t current_sum =  point_t::Zero(Dim);
        // recomputing waypoints q_i from derivative waypoints p_i. q_0 is the given constant.
        // then q_i = (sum( j = 0 -> j = i-1) p_j) /n+1
        n_wp.push_back(current_sum);
        for(typename t_point_t::const_iterator pit =  pts_.begin(); pit != pts_.end(); ++pit)
        {
            current_sum += *pit;
            n_wp.push_back(current_sum / new_degree);
        }
        bezier_curve_t integ(n_wp.begin(), n_wp.end(),T_min_, T_max_, mult_T_*(T_max_-T_min_));
        return integ.compute_primitive(order-1);
    }

    ///  \brief Evaluate the derivative order N of curve at time t.
    ///  If derivative is to be evaluated several times, it is
    ///  rather recommended to compute derived curve using compute_derivate.
    ///  \param order : order of derivative.
    ///  \param t : time when to evaluate the curve.
    ///  \return \f$\frac{d^Nx(t)}{dt^N}\f$ point corresponding on derived curve of order N at time t.
    ///
    virtual point_t derivate(const time_t t, const std::size_t order) const
    {
        bezier_curve_t deriv =compute_derivate(order);
        return deriv(t);
    }

    /// \brief Evaluate all Bernstein polynomes for a certain degree.
    /// A bezier curve with N control points is represented by : \f$x(t) = \sum_{i=0}^{N} B_i^N(t) P_i\f$
    /// with \f$ B_i^N(t) = \binom{N}{i}t^i (1-t)^{N-i} \f$.<br/>
    /// Warning: the horner scheme is about 100 times faster than this method.<br>
    /// This method will probably be removed in the future as the computation of bernstein polynomial is very costly.
    /// \param t : time when to evaluate the curve.
    /// \return \f$x(t)\f$ point corresponding on curve at time t.
    ///
    point_t evalBernstein(const Numeric t) const
    {
        const Numeric u = (t-T_min_)/(T_max_-T_min_);
        point_t res = point_t::Zero(Dim);
        typename t_point_t::const_iterator pts_it = pts_.begin();
        for(typename std::vector<Bern<Numeric> >::const_iterator cit = bernstein_.begin();
            cit !=bernstein_.end(); ++cit, ++pts_it)
        {
            res += cit->operator()(u) * (*pts_it);
        }
        return res*mult_T_;
    }

    /// \brief Evaluate all Bernstein polynomes for a certain degree using Horner's scheme.
    /// A bezier curve with N control points is expressed as : \f$x(t) = \sum_{i=0}^{N} B_i^N(t) P_i\f$.<br>
    /// To evaluate the position on curve at time t,we can apply the Horner's scheme : <br>
    /// \f$ x(t) = (1-t)^N(\sum_{i=0}^{N} \binom{N}{i} \frac{1-t}{t}^i P_i) \f$.<br>
    /// Horner's scheme : for a polynom of degree N expressed by : <br>
    /// \f$x(t) = a_0 + a_1t + a_2t^2 + ... + a_nt^n\f$
    /// where \f$number of additions = N\f$ / f$number of multiplication = N!\f$<br>
    /// Using Horner's method, the polynom is transformed into : <br>
    /// \f$x(t) = a_0 + t(a_1 + t(a_2+t(...))\f$ with N additions and multiplications.
    /// \param t : time when to evaluate the curve.
    /// \return \f$x(t)\f$ point corresponding on curve at time t.
    ///
    point_t evalHorner(const Numeric t) const
    {
        const Numeric u = (t-T_min_)/(T_max_-T_min_);
        typename t_point_t::const_iterator pts_it = pts_.begin();
        Numeric u_op, bc, tn;
        u_op = 1.0 - u;
        bc = 1;
        tn = 1;
        point_t tmp =(*pts_it)*u_op; ++pts_it;
        for(unsigned int i=1; i<degree_; i++, ++pts_it)
        {
            tn = tn*u;
            bc = bc*((num_t)(degree_-i+1))/i;
            tmp = (tmp + tn*bc*(*pts_it))*u_op;
        }
        return (tmp + tn*u*(*pts_it))*mult_T_;
    }

    const t_point_t& waypoints() const {return pts_;}

    /// \brief Evaluate the curve value at time t using deCasteljau algorithm.
    /// The algorithm will compute the \f$N-1\f$ centroids of parameters \f${t,1-t}\f$ of consecutive \f$N\f$ control points 
    /// of bezier curve, and perform it iteratively until getting one point in the list which will be the evaluation of bezier
    /// curve at time \f$t\f$.
    /// \param t : time when to evaluate the curve.
    /// \return \f$x(t)\f$ point corresponding on curve at time t.
    ///
    point_t evalDeCasteljau(const Numeric t) const {
        // normalize time :
        const Numeric u = (t-T_min_)/(T_max_-T_min_);
        t_point_t pts = deCasteljauReduction(waypoints(),u);
        while(pts.size() > 1)
        {
            pts = deCasteljauReduction(pts,u);
        }
        return pts[0]*mult_T_;
    }


    t_point_t deCasteljauReduction(const Numeric t) const{
        const Numeric u = (t-T_min_)/(T_max_-T_min_);
        return deCasteljauReduction(waypoints(),u);
    }

    /// \brief Compute de Casteljau's reduction of the given list of points at time t.
    /// For the list \f$pts\f$ of N points, compute a new list of points of size N-1 :<br>
    /// \f$<br>( pts[0]*(1-t)+pts[1], pts[1]*(1-t)+pts[2], ..., pts[0]*(N-2)+pts[N-1] )\f$<br>
    /// with t the time when to evaluate bezier curve.<br>\
    /// The new list contains centroid of parameters \f${t,1-t}\f$ of consecutive points in the list.
    /// \param pts : list of points.
    /// \param u   : NORMALIZED time when to evaluate the curve.
    /// \return reduced list of point (size of pts - 1).
    ///
    t_point_t deCasteljauReduction(const t_point_t& pts, const Numeric u) const{
        if(u < 0 || u > 1)
        {
            throw std::out_of_range("In deCasteljau reduction : u is not in [0;1]");
        }
        if(pts.size() == 1)
        {
            return pts;
        }

        t_point_t new_pts;
        for(cit_point_t cit = pts.begin() ; cit != (pts.end() - 1) ; ++cit)
        {
            new_pts.push_back((1-u) * (*cit) + u*(*(cit+1)));
        }
        return new_pts;
    }

    /// \brief Split the bezier curve in 2 at time t.
    /// \param t : list of points.
    /// \param u : unNormalized time.
    /// \return pair containing the first element of both bezier curve obtained.
    ///
    std::pair<bezier_curve_t,bezier_curve_t> split(const Numeric t){
        if (t == T_max_)
        {
            throw std::runtime_error("can't split curve, interval range is equal to original curve");
        }
        t_point_t wps_first(size_),wps_second(size_);
        const Numeric u = (t-T_min_)/(T_max_-T_min_);
        std::cout<<T_min_<<" and "<<T_max_<<" t="<<t<<" u="<<u<<std::endl;
        wps_first[0] = pts_.front();
        wps_second[degree_] = pts_.back();
        t_point_t casteljau_pts = waypoints();
        size_t id = 1;
        while(casteljau_pts.size() > 1)
        {
            casteljau_pts = deCasteljauReduction(casteljau_pts,u);
            wps_first[id] = casteljau_pts.front();
            wps_second[degree_-id] = casteljau_pts.back();
            ++id;
        }

        bezier_curve_t c_first(wps_first.begin(), wps_first.end(),T_min_,t,mult_T_);
        bezier_curve_t c_second(wps_second.begin(), wps_second.end(),t, T_max_,mult_T_);
        return std::make_pair(c_first,c_second);
    }

    bezier_curve_t extract(const Numeric t1, const Numeric t2){
        if(t1 < T_min_ || t1 > T_max_ || t2 < T_min_ || t2 > T_max_)
        {
            throw std::out_of_range("In Extract curve : times out of bounds");
        }
        if(t1 == T_min_ &&  t2 == T_max_)
        {
            return bezier_curve_t(waypoints().begin(), waypoints().end(), T_min_, T_max_, mult_T_);
        }
        if(t1 == T_min_)
        {
            return split(t2).first;
        }
        if(t2 == T_max_)
        {
            return split(t1).second;
        }

        std::pair<bezier_curve_t,bezier_curve_t> c_split = this->split(t1);
        return c_split.second.split(t2-t1).first;
    }

    private:
    template<typename In>
    t_point_t add_constraints(In PointsBegin, In PointsEnd, const curve_constraints_t& constraints)
    {
        t_point_t res;
        point_t P0, P1, P2, P_n_2, P_n_1, PN;
        P0 = *PointsBegin; PN = *(PointsEnd-1);
        P1    = P0+ constraints.init_vel / (num_t)degree_;
        P_n_1 = PN -constraints.end_vel  / (num_t)degree_;
        P2    = constraints.init_acc / (num_t)(degree_ * (degree_-1)) + 2* P1    - P0;
        P_n_2 = constraints.end_acc  / (num_t)(degree_ * (degree_-1)) + 2* P_n_1 - PN;

        res.push_back(P0);
        res.push_back(P1);
        res.push_back(P2);

        for(In it = PointsBegin+1; it != PointsEnd-1; ++it)
        {
            res.push_back(*it);
        }

        res.push_back(P_n_2);
        res.push_back(P_n_1);
        res.push_back(PN);
        return res;
    }

/*Operations*/

/*Helpers*/
    public:
    /// \brief Get the minimum time for which the curve is defined
    /// \return \f$t_{min}\f$, lower bound of time range.
    virtual time_t min() const{return T_min_;}
    /// \brief Get the maximum time for which the curve is defined.
    /// \return \f$t_{max}\f$, upper bound of time range.
    virtual time_t max() const{return T_max_;}
/*Helpers*/

	public:
    /// Starting time of cubic hermite spline : T_min_ is equal to first time of control points.
    /*const*/ time_t T_min_;
    /// Ending time of cubic hermite spline : T_max_ is equal to last time of control points.
    /*const*/ time_t T_max_;
    /*const*/ time_t mult_T_;
    /*const*/ std::size_t size_;
    /*const*/ std::size_t degree_;
    /*const*/ std::vector<Bern<Numeric> > bernstein_;
	
    private:
    t_point_t  pts_;

    public:
    static bezier_curve_t zero(const time_t T=1.)
    {
        std::vector<point_t> ts;
        ts.push_back(point_t::Zero(Dim));
        return bezier_curve_t(ts.begin(), ts.end(),0.,T);
    }
};
} // namespace curve
#endif //_CLASS_BEZIERCURVE

