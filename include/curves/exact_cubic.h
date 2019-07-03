/**
* \file exact_cubic.h
* \brief class allowing to create an Exact cubic spline.
* \author Steve T.
* \version 0.1
* \date 06/17/2013
*
* This file contains definitions for the ExactCubic class.
* Given a set of waypoints (x_i*) and timestep (t_i), it provides the unique set of
* cubic splines fulfulling those 4 restrictions :
* - x_i(t_i) = x_i* ; this means that the curve passes trough each waypoint
* - x_i(t_i+1) = x_i+1* ;
* - its derivative is continous at t_i+1
* - its acceleration is continous at t_i+1
* more details in paper "Task-Space Trajectories via Cubic Spline Optimization"
* By J. Zico Kolter and Andrew Y.ng (ICRA 2009)
*/


#ifndef _CLASS_EXACTCUBIC
#define _CLASS_EXACTCUBIC

#include "curve_abc.h"
#include "cubic_spline.h"
#include "quintic_spline.h"
#include "curve_constraint.h"

#include "MathDefs.h"

#include <functional>
#include <vector>

namespace curves
{
/// \class ExactCubic.
/// \brief Represents a set of cubic splines defining a continuous function 
/// crossing each of the waypoint given in its initialization.
///
template<typename Time= double, typename Numeric=Time, std::size_t Dim=3, bool Safe=false
, typename Point= Eigen::Matrix<Numeric, Dim, 1>, typename T_Point =std::vector<Point,Eigen::aligned_allocator<Point> >
, typename SplineBase=polynomial<Time, Numeric, Dim, Safe, Point, T_Point> >
struct exact_cubic : public curve_abc<Time, Numeric, Dim, Safe, Point>
{
    typedef Point   point_t;
    typedef T_Point t_point_t;
    typedef Eigen::Matrix<Numeric, Eigen::Dynamic, Eigen::Dynamic> MatrixX;
    typedef Eigen::Matrix<Numeric, 3, 3> Matrix3;
    typedef Time    time_t;
    typedef Numeric num_t;
    typedef SplineBase spline_t;
    typedef typename std::vector<spline_t> t_spline_t;
    typedef typename t_spline_t::iterator it_spline_t;
    typedef typename t_spline_t::const_iterator cit_spline_t;
    typedef curve_abc<Time, Numeric, Dim, Safe, Point> curve_abc_t;
    typedef curve_constraints<point_t> spline_constraints;

    /* Constructors - destructors */
    public:
    /// \brief Constructor.
    /// \param wayPointsBegin : an iterator pointing to the first element of a waypoint container.
    /// \param wayPointsEns   : an iterator pointing to the last element of a waypoint container.
    ///
    template<typename In>
    exact_cubic(In wayPointsBegin, In wayPointsEnd)
        : curve_abc_t(), subSplines_(computeWayPoints<In>(wayPointsBegin, wayPointsEnd)) {}

    /// \brief Constructor.
    /// \param wayPointsBegin : an iterator pointing to the first element of a waypoint container.
    /// \param wayPointsEns   : an iterator pointing to the last element of a waypoint container.
    /// \param constraints    : constraints on the init and end velocity / accelerations of the spline.
    ///
    template<typename In>
    exact_cubic(In wayPointsBegin, In wayPointsEnd, const spline_constraints& constraints)
        : curve_abc_t(), subSplines_(computeWayPoints<In>(wayPointsBegin, wayPointsEnd, constraints)) {}

    /// \brief Constructor.
    /// \param subSplines: vector of subsplines.
    exact_cubic(const t_spline_t& subSplines)
        : curve_abc_t(), subSplines_(subSplines) {}

    /// \brief Copy Constructor.
    exact_cubic(const exact_cubic& other)
        : curve_abc_t(), subSplines_(other.subSplines_) {}

    /// \brief Destructor.
    virtual ~exact_cubic(){}

    std::size_t getNumberSplines()
    {
        return subSplines_.size();
    }

    spline_t getSplineAt(std::size_t index)
    {
        return subSplines_.at(index);
    }

    private:
    /// \brief Compute polynom of exact cubic spline from waypoints.
    /// Compute the coefficients of polynom as in paper : "Task-Space Trajectories via Cubic Spline Optimization".<br>
    /// \f$x_i(t)=a_i+b_i(t-t_i)+c_i(t-t_i)^2\f$<br>
    /// with \f$a=x\f$, \f$H_1b=H_2x\f$, \f$c=H_3x+H_4b\f$, \f$d=H_5x+H_6b\f$.<br>
    /// The matrices \f$H\f$ are defined as in the paper in Appendix A.
    ///
    template<typename In>
    t_spline_t computeWayPoints(In wayPointsBegin, In wayPointsEnd) const
    {
        std::size_t const size(std::distance(wayPointsBegin, wayPointsEnd));
        if(Safe && size < 1)
        {
            throw std::length_error("size of waypoints must be superior to 0") ; // TODO
        }
        t_spline_t subSplines; subSplines.reserve(size);

        // refer to the paper to understand all this.
        MatrixX h1 = MatrixX::Zero(size, size);
        MatrixX h2 = MatrixX::Zero(size, size);
        MatrixX h3 = MatrixX::Zero(size, size);
        MatrixX h4 = MatrixX::Zero(size, size);
        MatrixX h5 = MatrixX::Zero(size, size);
        MatrixX h6 = MatrixX::Zero(size, size);

        MatrixX a =  MatrixX::Zero(size, Dim);
        MatrixX b =  MatrixX::Zero(size, Dim);
        MatrixX c =  MatrixX::Zero(size, Dim);
        MatrixX d =  MatrixX::Zero(size, Dim);
        MatrixX x =  MatrixX::Zero(size, Dim);


        In it(wayPointsBegin), next(wayPointsBegin);
        ++next;

        // Fill the matrices H as specified in the paper.
        for(std::size_t i(0); next != wayPointsEnd; ++next, ++it, ++i)
        {
            num_t const dTi((*next).first  - (*it).first);
            num_t const dTi_sqr(dTi * dTi);
            num_t const dTi_cube(dTi_sqr * dTi);
            // filling matrices values
            h3(i,i)   = -3 / dTi_sqr;
            h3(i,i+1) =  3 / dTi_sqr;
            h4(i,i)   = -2 / dTi;
            h4(i,i+1) = -1 / dTi;
            h5(i,i)   =  2 / dTi_cube;
            h5(i,i+1) = -2 / dTi_cube;
            h6(i,i)   =  1 / dTi_sqr;
            h6(i,i+1) =  1 / dTi_sqr;
            if( i+2 < size)
            {
                In it2(next); ++ it2;
                num_t const dTi_1((*it2).first - (*next).first);
                num_t const dTi_1sqr(dTi_1 * dTi_1);
                // this can be optimized but let's focus on clarity as long as not needed
                h1(i+1, i)   =  2 / dTi;
                h1(i+1, i+1) =  4 / dTi + 4 / dTi_1;
                h1(i+1, i+2) =  2 / dTi_1;
                h2(i+1, i)   = -6 / dTi_sqr;
                h2(i+1, i+1) = (6 / dTi_1sqr) - (6 / dTi_sqr);
                h2(i+1, i+2) =  6 / dTi_1sqr;
            }
            x.row(i)= (*it).second.transpose();
        }
        // adding last x
        x.row(size-1)= (*it).second.transpose();
        // Compute coefficients of polynom.
        a = x;
        PseudoInverse(h1);
        b = h1 * h2 * x; //h1 * b = h2 * x => b = (h1)^-1 * h2 * x
        c = h3 * x + h4 * b;
        d = h5 * x + h6 * b;
        // create splines along waypoints.
        it= wayPointsBegin, next=wayPointsBegin; ++ next;
        for(int i=0; next != wayPointsEnd; ++i, ++it, ++next)
        {
            subSplines.push_back(
                create_cubic<Time,Numeric,Dim,Safe,Point,T_Point>(a.row(i), b.row(i), c.row(i), d.row(i),(*it).first, (*next).first));
        }
        /*
        subSplines.push_back(
                create_cubic<Time,Numeric,Dim,Safe,Point,T_Point>(a.row(size-1), b.row(size-1), c.row(size-1), d.row(size-1), (*it).first, (*it).first));
        */
        return subSplines;
    }

    template<typename In>
    t_spline_t computeWayPoints(In wayPointsBegin, In wayPointsEnd, const spline_constraints& constraints) const
    {
        std::size_t const size(std::distance(wayPointsBegin, wayPointsEnd));
        if(Safe && size < 1) 
        {
            throw std::length_error("number of waypoints should be superior to one"); // TODO
        }
        t_spline_t subSplines; 
        subSplines.reserve(size-1);
        spline_constraints cons = constraints;
        In it(wayPointsBegin), next(wayPointsBegin), end(wayPointsEnd-1);
        ++next;
        for(std::size_t i(0); next != end; ++next, ++it, ++i)
        {
            compute_one_spline<In>(it, next, cons, subSplines);
        }
        compute_end_spline<In>(it, next,cons, subSplines);
        return subSplines;
    }

    template<typename In>
    void compute_one_spline(In wayPointsBegin, In wayPointsNext, spline_constraints& constraints, t_spline_t& subSplines) const
    {
        const point_t& a0 = wayPointsBegin->second, a1 = wayPointsNext->second;
        const point_t& b0 = constraints.init_vel , c0 = constraints.init_acc / 2.;
        const num_t& init_t = wayPointsBegin->first, end_t = wayPointsNext->first;
        const num_t dt = end_t - init_t, dt_2 = dt * dt, dt_3 = dt_2 * dt;
        const point_t d0 = (a1 - a0 - b0 * dt - c0 * dt_2) / dt_3;
        subSplines.push_back(create_cubic<Time,Numeric,Dim,Safe,Point,T_Point>
                             (a0,b0,c0,d0,init_t, end_t));
        constraints.init_vel = subSplines.back().derivate(end_t, 1);
        constraints.init_acc = subSplines.back().derivate(end_t, 2);
    }

    template<typename In>
    void compute_end_spline(In wayPointsBegin, In wayPointsNext, spline_constraints& constraints, t_spline_t& subSplines) const
    {
        const point_t& a0 = wayPointsBegin->second, a1 = wayPointsNext->second;
        const point_t& b0 = constraints.init_vel, b1 = constraints.end_vel,
                       c0 = constraints.init_acc / 2., c1 = constraints.end_acc;
        const num_t& init_t = wayPointsBegin->first, end_t = wayPointsNext->first;
        const num_t dt = end_t - init_t, dt_2 = dt * dt, dt_3 = dt_2 * dt, dt_4 = dt_3 * dt, dt_5 = dt_4 * dt;
        //solving a system of four linear eq with 4 unknows: d0, e0
        const point_t alpha_0 = a1 - a0 - b0 *dt -    c0 * dt_2;
        const point_t alpha_1 = b1 -      b0     - 2 *c0 * dt;
        const point_t alpha_2 = c1 -               2 *c0;
        const num_t x_d_0 = dt_3, x_d_1 = 3*dt_2, x_d_2 = 6*dt;
        const num_t x_e_0 = dt_4, x_e_1 = 4*dt_3, x_e_2 = 12*dt_2;
        const num_t x_f_0 = dt_5, x_f_1 = 5*dt_4, x_f_2 = 20*dt_3;

        point_t d, e, f;
        MatrixX rhs = MatrixX::Zero(3,Dim);
        rhs.row(0) = alpha_0; rhs.row(1) = alpha_1; rhs.row(2) = alpha_2;
        Matrix3 eq  = Matrix3::Zero(3,3);
        eq(0,0) = x_d_0; eq(0,1) = x_e_0; eq(0,2) = x_f_0;
        eq(1,0) = x_d_1; eq(1,1) = x_e_1; eq(1,2) = x_f_1;
        eq(2,0) = x_d_2; eq(2,1) = x_e_2; eq(2,2) = x_f_2;
        rhs = eq.inverse().eval() * rhs;
        d = rhs.row(0); e = rhs.row(1); f = rhs.row(2);

        subSplines.push_back(create_quintic<Time,Numeric,Dim,Safe,Point,T_Point>
                             (a0,b0,c0,d,e,f, init_t, end_t));
    }

    private:
    //exact_cubic& operator=(const exact_cubic&);
    /* Constructors - destructors */

    /*Operations*/
    public:
    ///  \brief Evaluation of the cubic spline at time t.
    ///  \param t : time when to evaluate the spline
    ///  \return \f$x(t)\f$ point corresponding on spline at time t.
    ///
    virtual point_t operator()(const time_t t) const
    {
        if(Safe && (t < subSplines_.front().min() || t > subSplines_.back().max()))
        {
            throw std::out_of_range("time t to evaluate should be in range [Tmin, Tmax] of the spline");
        }
        for(cit_spline_t it = subSplines_.begin(); it != subSplines_.end(); ++ it)
        {
            if( (t >= (it->min()) && t <= (it->max())) || it+1 == subSplines_.end())
            {
                return it->operator()(t);
            }
        }
        // this should not happen
        throw std::runtime_error("Exact cubic evaluation failed; t is outside bounds");
    }

    ///  \brief Evaluate the derivative of order N of spline at time t.
    ///  \param t : time when to evaluate the spline.
    ///  \param order : order of derivative.
    ///  \return \f$\frac{d^Nx(t)}{dt^N}\f$ point corresponding on derivative spline of order N at time t.
    ///
    virtual point_t derivate(const time_t t, const std::size_t order) const
    {
        if(Safe && (t < subSplines_.front().min() || t > subSplines_.back().max()))
        {
            throw std::out_of_range("time t to evaluate should be in range [Tmin, Tmax] of the spline");
        }
        for(cit_spline_t it = subSplines_.begin(); it != subSplines_.end(); ++ it)
        {
            if( (t >= (it->min()) && t <= (it->max())) || it+1 == subSplines_.end())
            {
                return it->derivate(t, order);
            }
        }
        // this should not happen
        throw std::runtime_error("Exact cubic evaluation failed; t is outside bounds");
    }
    /*Operations*/

    /*Helpers*/
    public:
    /// \brief Get the minimum time for which the curve is defined
    /// \return \f$t_{min}\f$ lower bound of time range.
    num_t virtual min() const{return subSplines_.front().min();}
    /// \brief Get the maximum time for which the curve is defined.
    /// \return \f$t_{max}\f$ upper bound of time range.
    num_t virtual max() const{return subSplines_.back().max();}
    /*Helpers*/

    /*Attributes*/
    public:
    t_spline_t subSplines_; // const
    /*Attributes*/
};
} // namespace curves
#endif //_CLASS_EXACTCUBIC

