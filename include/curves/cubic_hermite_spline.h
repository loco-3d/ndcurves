#ifndef _CLASS_CUBICHERMITESPLINE
#define _CLASS_CUBICHERMITESPLINE

#include "curve_abc.h"
#include "bernstein.h"
#include "curve_constraint.h"

#include "MathDefs.h"

#include <vector>
#include <stdexcept>

#include <iostream>

namespace curves
{
/// \class CubicHermiteSpline.
/// \brief Represents a set of cubic hermite splines defining a continuous function \f$x(t)\f$.
/// A hermite cubic spline is a minimal degree polynom interpolating a function in two 
/// points \f$P_i\f$ and \f$P_{i+1}\f$ with its tangent \f$m_i\f$ and \f$m_{i+1}\f$.
/// A hermite cubic spline :
/// - crosses each of the waypoint given in its initialization (\f$P_0\f$, \f$P_1\f$,...,\f$P_N\f$).
/// - has its derivate on \f$P_i\f$ and \f$P_{i+1}\f$ are \f$x'(t_{P_i}) = m_i\f$ and \f$x'(t_{P_{i+1}}) = m_{i+1}\f$.
///
template<typename Time= double, typename Numeric=Time, std::size_t Dim=3, bool Safe=false
, typename Point= Eigen::Matrix<Numeric, Dim, 1>
, typename Tangent= Eigen::Matrix<Numeric, Dim, 1>
, typename Pair_point_tangent= std::pair<Point, Tangent>
, typename Vector_pair= std::vector< Pair_point_tangent ,Eigen::aligned_allocator<Point> > 
>
struct cubic_hermite_spline : public curve_abc<Time, Numeric, Dim, Safe, Point>
{

    typedef int Index;
    typedef std::vector<Time> Vector_time;

    /*Attributes*/
    public:
    Vector_pair control_points_; // Vector of pair < Point, Tangent > , Dim = Size_.
    Vector_time time_control_points_; // Time corresponding to each control point, Dim = Size_.
    Vector_time duration_splines_; // its length should be duration_splines_[0] = t_PO_P1, duration_splines_[1] = t_P1_P2, Dim = Size_-1.
    /*const*/ Time T_min_;
    /*const*/ Time T_max_;
    /*const*/ std::size_t size_;
    /*Attributes*/
    
    public:
	/// \brief Constructor.
	/// \param wayPointsBegin : an iterator pointing to the first element of a waypoint container.
	/// \param wayPointsEns   : an iterator pointing to the last element of a waypoint container.
    ///
	template<typename In>
	cubic_hermite_spline(In PairsBegin, In PairsEnd)
	: T_min_(0.)
    , T_max_(1.)
    , size_(std::distance(PairsBegin, PairsEnd))
	{
		// Check size of pairs container.
        std::size_t const size(std::distance(PairsBegin, PairsEnd));
        if(Safe && size < 1)
        {
            throw std::length_error("can't create cubic_hermite_spline, number of pairs is inferior to 2.");
        }
        // Push all pairs in controlPoints
        In it(PairsBegin);
        for(; it != PairsEnd; ++it)
        {
            control_points_.push_back(*it);
        }
        setTimeSplinesDefault();
	}

	/// \brief Destructor.
    virtual ~cubic_hermite_spline(){}

    /*Operations*/
	public:
	virtual Point operator()(const Time t) const
    {
        if(Safe &! (T_min_ <= t && t <= T_max_))
        {
			throw std::out_of_range("can't evaluate bezier curve, out of range");
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
    ///  \return \f$\frac{d^Nx(t)}{dt^N}\f$, point corresponding on derivative spline of order N at time t.
    ///
    virtual Point derivate(const Time t, const std::size_t order) const
    {   
        return evalCubicHermiteSpline(t, order);
    }

    /// \brief Set duration of each spline.
    /// Set duration of each spline, exemple : time_control_points[0] = duration of first spline, 
    /// time_control_points[0] = duration of second spline (...).
    /// \param time_control_points : Vector containing duration of each spline.
    void setTimeSplines(const Vector_time & time_control_points)
	{
        time_control_points_ = time_control_points;
        T_min_ = time_control_points_.front();
        T_max_ = time_control_points_.back();
		assert(time_control_points.size() == size());
        compute_duration_splines();
        if (!check_duration_splines())
        {
            throw std::logic_error("time_splines not monotonous, all spline duration should be superior to 0");
        }
	}

    /// \brief Set duration by default of each spline.
    /// Set duration by default (duration = T / Number of control points) for each spline.
    ///
    void setTimeSplinesDefault()
    {
        time_control_points_.clear();
        T_min_ = 0.;
        T_max_ = 1.;
        Time timestep = (T_max_- T_min_) / (control_points_.size()-1);
        Time time = 0.;
        Index i = 0;
        for (i=0; i<size(); i++)
        {
            //time_control_points_.push_back(time);
            time_control_points_.push_back(time);
            time += timestep;
        }
        compute_duration_splines();
    }

    /// \brief Returns the index of the interval for the interpolation.
    ///
    Index findInterval(const Numeric t) const
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

        Index left_id = 0;
        Index right_id = size_-1;
        while(left_id <= right_id)
        {
            const Index middle_id = left_id + (right_id - left_id)/2;
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

    /*
        Vector_pair control_points_; // Vector of pair < Point, Tangent >
        Vector_time time_control_points_; // Time corresponding to each control point.
        Vector_time duration_splines_;
    */
    // See how to evaluate cubic hermite spline
    Point evalCubicHermiteSpline(const Numeric t, std::size_t order_derivated) const
    {
        // TO DO
        const Index id = findInterval(t);

        // ID is on the last control point
        if(id == size_-1)
        {
            if (order_derivated==0)
            {
                return control_points_.back().first;
            }
            else if (order_derivated==1)
            {
                return control_points_.back().second;
            }
            else
            {
                return control_points_.back().first*0.; // To modify, create a new Tangent ininitialized with 0.
            }
        }

        const Pair_point_tangent Pair0 = control_points_.at(id);
        const Pair_point_tangent Pair1 = control_points_.at(id+1);

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

        assert(0. <= alpha <= 1. && "alpha must be in [0,1]");

        Numeric h00, h10, h01, h11;
        evalCoeffs(alpha,h00,h10,h01,h11,order_derivated);
        //std::cout << "for val t="<<t<<" coef : h00="<<h00<<" h10="<<h10<<" h01="<<h01<<" h11="<<h11<<std::endl;

        Point p_ = (h00 * Pair0.first + h10 * Pair0.second + h01 * Pair1.first + h11 * Pair1.second);
        return p_;
    }

    
    /// \brief Returns the number of points contained in the trajectory.
    ///
    Index size() const { return size_; }

    /// \brief Returns the number of intervals (splines) contained in the trajectory.
    ///
    Index numIntervals() const { return size()-1; }


	private:

    /// \brief compute duration of each spline, ex : \f$t_{P_0\_P_1}=t_{P_1}-t_{P_0}f$.
    void compute_duration_splines() {
        duration_splines_.clear();
        Time actual_time;
        Time prev_time = *(time_control_points_.begin());
        Index i = 0;
        for (i=0; i<size()-1; i++)
        {
            actual_time = time_control_points_.at(i+1);
            duration_splines_.push_back(actual_time-prev_time);
            prev_time = actual_time;
        }
    }

    /// \brief Check if the absicca
    bool check_duration_splines() const
    {
        Index i = 0;
        bool is_positive = true;
        while (is_positive && i<duration_splines_.size())
        {
            is_positive = (duration_splines_.at(i)>0.);
            i++;
        }
        return is_positive;
    }

    static void evalCoeffs(const Numeric t, Numeric & h00, Numeric & h10, Numeric & h01, Numeric & h11, std::size_t order_derivated)
    {
        Numeric t_square = t*t;
        Numeric t_cube = t_square*t;
        if (order_derivated==0)
        {
            h00 =  2*t_cube     -3*t_square     +1.;
            h10 =  t_cube       -2*t_square     +t;
            h01 = -2*t_cube     +3*t_square;
            h11 =  t_cube       -t_square;
        }
        else if (order_derivated==1)
        {
            h00 =  6*t_square   -6*t;
            h10 =  3*t_square   -4*t    +1.;
            h01 = -6*t_square   +6*t;
            h11 =  3*t_square   -2*t;
        }
        else if (order_derivated==2)
        {
            h00 =  12*t     -6.;
            h10 =  6*t      -4.;  
            h01 = -12*t     +6.;
            h11 =  6*t      -2.;
        }
        else if (order_derivated==3)
        {
            h00 =  12.;
            h10 =  6.;  
            h01 = -12.;
            h11 =  0.;
        }
        else 
        {
            h00 =  0.;
            h10 =  0.;  
            h01 =  0.;
            h11 =  0.;
        }
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

};

} // namespace curve
#endif //_CLASS_CUBICHERMITESPLINE