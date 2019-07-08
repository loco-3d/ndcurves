/**
* \file cubic_hermite_spline.h
* \brief class allowing to create a cubic hermite spline of any dimension.
* \author Justin Carpentier <jcarpent@laas.fr> modified by Jason Chemin <jchemin@laas.fr>
* \date 05/2019
*/

#ifndef _CLASS_CUBICHERMITESPLINE
#define _CLASS_CUBICHERMITESPLINE

#include "curve_abc.h"
#include "curve_constraint.h"

#include "MathDefs.h"

#include <vector>
#include <stdexcept>

#include <iostream>

namespace curves
{
/// \class CubicHermiteSpline.
/// \brief Represents a set of cubic hermite splines defining a continuous function \f$p(t)\f$.
/// A hermite cubic spline is a minimal degree polynom interpolating a function in two 
/// points \f$P_i\f$ and \f$P_{i+1}\f$ with its tangent \f$m_i\f$ and \f$m_{i+1}\f$.<br>
/// A hermite cubic spline :
/// - crosses each of the waypoint given in its initialization (\f$P_0\f$, \f$P_1\f$,...,\f$P_N\f$).
/// - has its derivatives on \f$P_i\f$ and \f$P_{i+1}\f$ are \f$p'(t_{P_i}) = m_i\f$ and \f$p'(t_{P_{i+1}}) = m_{i+1}\f$.
///
template<typename Time= double, typename Numeric=Time, std::size_t Dim=3, bool Safe=false
, typename Point= Eigen::Matrix<Numeric, Eigen::Dynamic, 1>
, typename Tangent= Eigen::Matrix<Numeric, Eigen::Dynamic, 1>
>
struct cubic_hermite_spline : public curve_abc<Time, Numeric, Safe, Point>
{
    typedef std::pair<Point, Tangent> pair_point_tangent_t; 
    typedef std::vector< pair_point_tangent_t ,Eigen::aligned_allocator<Point> > t_pair_point_tangent_t;
    typedef std::vector<Time> vector_time_t;
    typedef Numeric num_t;
    
    public:
    /// \brief Constructor.
    /// \param wayPointsBegin : an iterator pointing to the first element of a pair(position, derivative) container.
    /// \param wayPointsEns   : an iterator pointing to the last  element of a pair(position, derivative) container.
    /// \param time_control_points : vector containing time for each waypoint.
    ///
    template<typename In>
    cubic_hermite_spline(In PairsBegin, In PairsEnd, const vector_time_t & time_control_points)
    : degree_(3)
    {
        // Check size of pairs container.
        std::size_t const size(std::distance(PairsBegin, PairsEnd));
        size_ = size;
        if(Safe && size < 1)
        {
            throw std::length_error("can not create cubic_hermite_spline, number of pairs is inferior to 2.");
        }
        // Push all pairs in controlPoints
        In it(PairsBegin);
        for(; it != PairsEnd; ++it)
        {
            control_points_.push_back(*it);
        }
        setTime(time_control_points);
    }

    cubic_hermite_spline(const cubic_hermite_spline& other)
        : t_min_(other.t_min_), t_max_(other.t_max_), size_(other.size_), degree_(other.degree_), 
          control_points_(other.control_points_), time_control_points_(other.time_control_points_),
          duration_splines_(other.duration_splines_)
          {}

    /// \brief Destructor.
    virtual ~cubic_hermite_spline(){}

    /*Operations*/
    public:
    ///  \brief Evaluation of the cubic hermite spline at time t.
    ///  \param t : time when to evaluate the spline.
    ///  \return \f$p(t)\f$ point corresponding on spline at time t.
    ///
    virtual Point operator()(const Time t) const
    {
        if(Safe &! (t_min_ <= t && t <= t_max_))
        {
            throw std::invalid_argument("can't evaluate cubic hermite spline, out of range");
        }
        if (size_ == 1)
        {
            return control_points_.front().first;
        }
        else
        {
            return evalCubicHermiteSpline(t, 0);
        }
    }

    ///  \brief Evaluate the derivative of order N of spline at time t.
    ///  \param t : time when to evaluate the spline.
    ///  \param order : order of derivative.
    ///  \return \f$\frac{d^Np(t)}{dt^N}\f$ point corresponding on derivative spline of order N at time t.
    ///
    virtual Point derivate(const Time t, const std::size_t order) const
    {   
        return evalCubicHermiteSpline(t, order);
    }

    /// \brief Set time of each control point of cubic hermite spline.
    /// Set duration of each spline, Exemple : \f$( 0., 0.5, 0.9, ..., 4.5 )\f$ with 
    /// values corresponding to times for \f$P_0, P_1, P_2, ..., P_N\f$ respectively.<br>
    /// \param time_control_points : Vector containing time for each control point.
    ///
    void setTime(const vector_time_t & time_control_points)
    {
        time_control_points_ = time_control_points;
        t_min_ = time_control_points_.front();
        t_max_ = time_control_points_.back();
        if (time_control_points.size() != size())
        {
            throw std::length_error("size of time control points should be equal to number of control points");
        }
        computeDurationSplines();
        if (!checkDurationSplines())
        {
            throw std::invalid_argument("time_splines not monotonous, all spline duration should be superior to 0");
        }
    }

    /// \brief Get vector of pair (positition, derivative) corresponding to control points.
    /// \return vector containing control points.
    ///
    t_pair_point_tangent_t getControlPoints()
    {
        return control_points_;
    }

    /// \brief Get vector of Time corresponding to Time for each control point.
    /// \return vector containing time of each control point.
    ///
    vector_time_t getTime()
    {
        return time_control_points_;
    }

    /// \brief Get number of control points contained in the trajectory.
    /// \return number of control points.
    ///
    std::size_t size() const { return size_; }

    /// \brief Get number of intervals (subsplines) contained in the trajectory.
    /// \return number of intervals (subsplines).
    ///
    std::size_t numIntervals() const { return size()-1; }


    /// \brief Evaluate value of cubic hermite spline or its derivate at specified order at time \f$t\f$.
    /// A cubic hermite spline on unit interval \f$[0,1]\f$ and given two control points defined by 
    /// their position and derivative \f$\{p_0,m_0\}\f$ and \f$\{p_1,m_1\}\f$, is defined by the polynom : <br>
    ///     \f$p(t)=h_{00}(t)P_0 + h_{10}(t)m_0 + h_{01}(t)p_1 + h_{11}(t)m_1\f$<br>
    /// To extend this formula to a cubic hermite spline on any arbitrary interval, 
    /// we define \f$\alpha=(t-t_0)/(t_1-t_0) \in [0, 1]\f$ where \f$t \in [t_0, t_1]\f$.<br>
    /// Polynom \f$p(t) \in [t_0, t_1]\f$ becomes \f$p(\alpha) \in [0, 1]\f$
    /// and \f$p(\alpha) = p((t-t_0)/(t_1-t_0))\f$.
    /// \param t : time when to evaluate the curve.
    /// \param degree_derivative : Order of derivate of cubic hermite spline (set value to 0 if you do not want derivate)
    /// \return point corresponding \f$p(t)\f$ on spline at time t or its derivate order N \f$\frac{d^Np(t)}{dt^N}\f$.
    ///
    Point evalCubicHermiteSpline(const Numeric t, std::size_t degree_derivative) const
    {
        const std::size_t id = findInterval(t);
        // ID is on the last control point
        if(id == size_-1)
        {
            if (degree_derivative==0)
            {
                return control_points_.back().first;
            }
            else if (degree_derivative==1)
            {
                return control_points_.back().second;
            }
            else
            {
                return control_points_.back().first*0.; // To modify, create a new Tangent ininitialized with 0.
            }
        }
        const pair_point_tangent_t Pair0 = control_points_.at(id);
        const pair_point_tangent_t Pair1 = control_points_.at(id+1);
        const Time & t0 = time_control_points_[id];
        const Time & t1 = time_control_points_[id+1]; 
        // Polynom for a cubic hermite spline defined on [0., 1.] is : 
        //      p(t) = h00(t)*p0 + h10(t)*m0 + h01(t)*p1 + h11(t)*m1 with t in [0., 1.]
        //
        // For a cubic hermite spline defined on [t0, t1], we define alpha=(t-t0)/(t1-t0) in [0., 1.].
        // Polynom p(t) defined on [t0, t1] becomes p(alpha) defined on [0., 1.]
        //      p(alpha) = p((t-t0)/(t1-t0))
        //
        const Time dt = (t1-t0);
        const Time alpha = (t - t0)/dt;
        assert(0. <= alpha && alpha <= 1. && "alpha must be in [0,1]");
        Numeric h00, h10, h01, h11;
        evalCoeffs(alpha,h00,h10,h01,h11,degree_derivative);
        //std::cout << "for val t="<<t<<" alpha="<<alpha<<" coef : h00="<<h00<<" h10="<<h10<<" h01="<<h01<<" h11="<<h11<<std::endl;
        Point p_ = (h00 * Pair0.first + h10 * dt * Pair0.second + h01 * Pair1.first + h11 * dt * Pair1.second);
        // if derivative, divide by dt^degree_derivative
        for (std::size_t i=0; i<degree_derivative; i++)
        {
            p_ /= dt;
        }
        return p_;
    }

    /// \brief Evaluate coefficient for polynom of cubic hermite spline.
    /// Coefficients of polynom :<br>
    ///  - \f$h00(t)=2t^3-3t^2+1\f$;
    ///  - \f$h10(t)=t^3-2t^2+t\f$;
    ///  - \f$h01(t)=-2t^3+3t^2\f$;
    ///  - \f$h11(t)=t^3-t^2\f$.<br>
    /// From it, we can calculate their derivate order N : 
    /// \f$\frac{d^Nh00(t)}{dt^N}\f$, \f$\frac{d^Nh10(t)}{dt^N}\f$,\f$\frac{d^Nh01(t)}{dt^N}\f$, \f$\frac{d^Nh11(t)}{dt^N}\f$.
    /// \param t : time to calculate coefficients.
    /// \param h00 : variable to store value of coefficient.
    /// \param h10 : variable to store value of coefficient.
    /// \param h01 : variable to store value of coefficient.
    /// \param h11 : variable to store value of coefficient.
    /// \param degree_derivative : order of derivative.
    ///
    static void evalCoeffs(const Numeric t, Numeric & h00, Numeric & h10, Numeric & h01, Numeric & h11, std::size_t degree_derivative)
    {
        Numeric t_square = t*t;
        Numeric t_cube = t_square*t;
        if (degree_derivative==0)
        {
            h00 =  2*t_cube     -3*t_square     +1.;
            h10 =  t_cube       -2*t_square     +t;
            h01 = -2*t_cube     +3*t_square;
            h11 =  t_cube       -t_square;
        }
        else if (degree_derivative==1)
        {
            h00 =  6*t_square   -6*t;
            h10 =  3*t_square   -4*t    +1.;
            h01 = -6*t_square   +6*t;
            h11 =  3*t_square   -2*t;
        }
        else if (degree_derivative==2)
        {
            h00 =  12*t     -6.;
            h10 =  6*t      -4.;  
            h01 = -12*t     +6.;
            h11 =  6*t      -2.;
        }
        else if (degree_derivative==3)
        {
            h00 =  12.;
            h10 =  6.;  
            h01 = -12.;
            h11 =  6.;
        }
        else 
        {
            h00 =  0.;
            h10 =  0.;  
            h01 =  0.;
            h11 =  0.;
        }
    }

    private:
    /// \brief Get index of the interval (subspline) corresponding to time t for the interpolation.
    /// \param t : time where to look for interval.
    /// \return Index of interval for time t.
    ///
    std::size_t findInterval(const Numeric t) const
    {
        // time before first control point time.
        if(t < time_control_points_[0])
        {
            return 0;
        }
        // time is after last control point time
        if(t > time_control_points_[size_-1])
        {
            return size_-1;
        }

        std::size_t left_id = 0;
        std::size_t right_id = size_-1;
        while(left_id <= right_id)
        {
            const std::size_t middle_id = left_id + (right_id - left_id)/2;
            if(time_control_points_.at(middle_id) < t)
            {
                left_id = middle_id+1;
            }
            else if(time_control_points_.at(middle_id) > t)
            {
                right_id = middle_id-1;
            }
            else
            {
                return middle_id;
            }
        }
        return left_id-1;
    }

    /// \brief compute duration of each spline.
    /// For N control points with time \f$T_{P_0}, T_{P_1}, T_{P_2}, ..., T_{P_N}\f$ respectively,
    /// Duration of each subspline is : ( T_{P_1}-T_{P_0}, T_{P_2}-T_{P_1}, ..., T_{P_N}-T_{P_{N-1} ).
    ///
    void computeDurationSplines() {
        duration_splines_.clear();
        Time actual_time;
        Time prev_time = *(time_control_points_.begin());
        std::size_t i = 0;
        for (i=0; i<size()-1; i++)
        {
            actual_time = time_control_points_.at(i+1);
            duration_splines_.push_back(actual_time-prev_time);
            prev_time = actual_time;
        }
    }

    /// \brief Check if duration of each subspline is strictly positive.
    /// \return true if all duration of strictly positive, false otherwise.
    ///
    bool checkDurationSplines() const
    {
        std::size_t i = 0;
        bool is_positive = true;
        while (is_positive && i<duration_splines_.size())
        {
            is_positive = (duration_splines_.at(i)>0.);
            i++;
        }
        return is_positive;
    }
    /*Operations*/  

    /*Helpers*/
    public:
    /// \brief Get the minimum time for which the curve is defined
    /// \return \f$t_{min}\f$, lower bound of time range.
    Time virtual min() const{return time_control_points_.front();}
    /// \brief Get the maximum time for which the curve is defined.
    /// \return \f$t_{max}\f$, upper bound of time range.
    Time virtual max() const{return time_control_points_.back();}
    /*Helpers*/

    /*Attributes*/
    /// Vector of pair < Point, Tangent >.
    t_pair_point_tangent_t control_points_;
    /// Vector of Time corresponding to time of each N control points : time at \f$P_0, P_1, P_2, ..., P_N\f$.
    /// Exemple : \f$( 0., 0.5, 0.9, ..., 4.5 )\f$ with values corresponding to times for \f$P_0, P_1, P_2, ..., P_N\f$ respectively.
    vector_time_t time_control_points_;

    /// Vector of Time corresponding to time duration of each subspline.<br>
    /// For N control points with time \f$T_{P_0}, T_{P_1}, T_{P_2}, ..., T_{P_N}\f$ respectively,
    /// duration of each subspline is : ( T_{P_1}-T_{P_0}, T_{P_2}-T_{P_1}, ..., T_{P_N}-T_{P_{N-1} )<br>
    /// It contains \f$N-1\f$ durations. 
    vector_time_t duration_splines_;
    /// Starting time of cubic hermite spline : t_min_ is equal to first time of control points.
    /*const*/ Time t_min_;
    /// Ending time of cubic hermite spline : t_max_ is equal to last time of control points.
    /*const*/ Time t_max_;
    /// Number of control points (pairs).
    std::size_t size_;
    /// Degree (Cubic so degree 3)
    std::size_t degree_;
    /*Attributes*/

};

} // namespace curve
#endif //_CLASS_CUBICHERMITESPLINE