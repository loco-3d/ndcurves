#ifndef _CLASS_PIECEWISE_CURVE
#define _CLASS_PIECEWISE_CURVE

#include "curve_abc.h"
#include "polynomial.h"
#include "curve_conversion.h"


namespace curves
{
/// \class PiecewiseCurve.
/// 
///
template<typename Time= double, typename Numeric=Time, std::size_t Dim=3, bool Safe=false,
     typename Point= Eigen::Matrix<Numeric, Dim, 1>, 
     typename T_Point= std::vector<Point,Eigen::aligned_allocator<Point> > 
     >
struct piecewise_polynomial_curve : public curve_abc<Time, Numeric, Dim, Safe, Point>
{

	typedef Point 	point_t;
	typedef T_Point t_point_t;
	typedef Time 	time_t;
    typedef Numeric	num_t;
    typedef int Index;

    //typedef polynomial  <double, double, 3, true, point_t, t_point_t> polynomial_t;
    typedef polynomial  <double, double, 3, true, point_t, t_point_t> polynomial_t;
    typedef typename std::vector < polynomial_t > t_polynomial_t;
    typedef std::vector<Time> t_vector_time_t;

	public:

	/// \brief Constructor.
    /// Initialize a piecewise polynomial curve by giving the first polynomial curve.
    /// \param pol   : a polynomial curve.
    ///
	piecewise_polynomial_curve(const polynomial_t& pol)
	{
		size_ = 0;
		T_min_ = pol.min();
		time_polynomial_curves_.push_back(T_min_);
		add_polynomial_curve(pol);
	}

	virtual ~piecewise_polynomial_curve(){}

	virtual Point operator()(const Time t) const
    {
        if(Safe &! (T_min_ <= t && t <= T_max_))
        {
			throw std::out_of_range("can't evaluate piecewise curve, out of range");
        }
        return polynomial_curves_.at(findInterval(t))(t);
    }

    ///  \brief Evaluate the derivative of order N of curve at time t.
    ///  \param t : time when to evaluate the spline.
    ///  \param order : order of derivative.
    ///  \return \f$\frac{d^Np(t)}{dt^N}\f$ point corresponding on derivative spline of order N at time t.
    ///
    virtual Point derivate(const Time t, const std::size_t order) const
    {   
        if(Safe &! (T_min_ <= t && t <= T_max_))
        {
			throw std::out_of_range("can't evaluate piecewise curve, out of range");
        }
        return (polynomial_curves_.at(findInterval(t))).derivate(t, order);
    }

	void add_polynomial_curve(polynomial_t pol)
	{
		// Check time continuity : Begin time of pol must be equal to T_max_ of actual piecewise curve.
		if (size_!=0 && pol.min()!=T_max_)
		{
			throw std::out_of_range("Can not add new Polynom to PiecewiseCurve : time discontinuity between T_max_ and pol.min()");
		}
		polynomial_curves_.push_back(pol);
		size_ = polynomial_curves_.size();
		T_max_ = pol.max();
		time_polynomial_curves_.push_back(T_max_);
	}

	bool isContinuous(const std::size_t order)
	{
		bool isContinuous =true;
    	Index i=0;
    	point_t value_end, value_start;
    	while (isContinuous && i<(size_-1))
    	{
    		polynomial_t current = polynomial_curves_.at(i);
    		polynomial_t next = polynomial_curves_.at(i+1);

    		if (order == 0)
    		{
    			value_end = current(current.max());
				value_start = next(next.min());
    		}
    		else
    		{
				value_end = current.derivate(current.max(),order);
				value_start = next.derivate(next.min(),order);
			}

    		if ((value_end-value_start).norm() > margin)
    		{
    			isContinuous = false;
    		}
    		i++;
    	}
    	return isContinuous;
	}

	private:

	/// \brief Get index of the interval corresponding to time t for the interpolation.
    /// \param t : time where to look for interval.
    /// \return Index of interval for time t.
    ///
    Index findInterval(const Numeric t) const
    {	
        // time before first control point time.
        if(t < time_polynomial_curves_[0])
        {
            return 0;
        }
        // time is after last control point time
        if(t > time_polynomial_curves_[size_-1])
        {
            return size_-1;
        }

        Index left_id = 0;
        Index right_id = size_-1;
        while(left_id <= right_id)
        {
            const Index middle_id = left_id + (right_id - left_id)/2;
            if(time_polynomial_curves_.at(middle_id) < t)
            {
                left_id = middle_id+1;
            }
            else if(time_polynomial_curves_.at(middle_id) > t)
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

    /*Helpers*/
	public:
    /// \brief Get the minimum time for which the curve is defined
    /// \return \f$t_{min}\f$, lower bound of time range.
    Time virtual min() const{return T_min_;}
    /// \brief Get the maximum time for which the curve is defined.
    /// \return \f$t_{max}\f$, upper bound of time range.
    Time virtual max() const{return T_max_;}
	/*Helpers*/

    /* Variables */
	t_polynomial_t polynomial_curves_; // for curves 0/1/2 : [ curve0, curve1, curve2 ]
	t_vector_time_t time_polynomial_curves_; // for curves 0/1/2 : [ Tmin0, Tmax0,Tmax1,Tmax2 ]
	Numeric size_; // Number of segments in piecewise curve
	Time T_min_, T_max_;
	const double margin = 0.001;
};

} // end namespace


#endif // _CLASS_PIECEWISE_CURVE