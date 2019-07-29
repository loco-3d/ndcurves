/**
* \file piecewise_curve.h
* \brief class allowing to create a piecewise curve.
* \author Jason C.
* \date 05/2019
*/

#ifndef _CLASS_PIECEWISE_CURVE
#define _CLASS_PIECEWISE_CURVE

#include "curve_abc.h"
#include "curve_conversion.h"

#include <typeinfo>

namespace curves
{
/// \class PiecewiseCurve.
/// \brief Represent a piecewise curve. We can add some new curve,
///        but the starting time of the curve to add should be equal to the ending time of the actual
///        piecewise_curve.<br>\ Example : A piecewise curve composed of three curves cf0, 
///        cf1 and cf2 where cf0 is defined between \f$[T0_{min},T0_{max}]\f$, cf1 between 
///        \f$[T0_{max},T1_{max}]\f$ and cf2 between \f$[T1_{max},T2_{max}]\f$.
///        On the piecewise polynomial curve, cf0 is located between \f$[T0_{min},T0_{max}[\f$,
///        cf1 between \f$[T0_{max},T1_{max}[\f$ and cf2 between \f$[T1_{max},T2_{max}]\f$.
///
template<typename Time= double, typename Numeric=Time, std::size_t Dim=3, bool Safe=false,
     typename Point= Eigen::Matrix<Numeric, Eigen::Dynamic, 1>, 
     typename T_Point= std::vector<Point,Eigen::aligned_allocator<Point> >,
     typename Curve= curve_abc<Time, Numeric, Safe, Point>
     >
struct piecewise_curve : public curve_abc<Time, Numeric, Safe, Point>
{
    typedef Point   point_t;
    typedef T_Point t_point_t;
    typedef Time    time_t;
    typedef Numeric num_t;
    typedef Curve curve_t;
    typedef typename std::vector < curve_t > t_curve_t;
    typedef typename std::vector< Time > t_time_t;

    public:

    /// \brief Constructor.
    /// Initialize a piecewise curve by giving the first curve.
    /// \param pol   : a polynomial curve.
    ///
    piecewise_curve(const curve_t& cf)
    {
        size_ = 0;
        add_curve(cf);
        time_curves_.push_back(cf.min());
        T_min_ = cf.min();
    }

    virtual ~piecewise_curve(){}

    virtual Point operator()(const Time t) const
    {
        if(Safe &! (T_min_ <= t && t <= T_max_))
        {
            //std::cout<<"[Min,Max]=["<<T_min_<<","<<T_max_<<"]"<<" t="<<t<<std::endl;
            throw std::out_of_range("can't evaluate piecewise curve, out of range");
        }
        return curves_.at(find_interval(t))(t);
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
            throw std::invalid_argument("can't evaluate piecewise curve, out of range");
        }
        return (curves_.at(find_interval(t))).derivate(t, order);
    }

    ///  \brief Add a new curve to piecewise curve, which should be defined in \f$[T_{min},T_{max}]\f$ where \f$T_{min}\f$
    ///         is equal to \f$T_{max}\f$ of the actual piecewise curve.
    ///  \param cf : curve to add.
    ///
    void add_curve(const curve_t& cf)
    {
        // Check time continuity : Beginning time of pol must be equal to T_max_ of actual piecewise curve.
        if (size_!=0 && !(fabs(cf.min()-T_max_)<MARGIN))
        {
            throw std::invalid_argument("Can not add new Polynom to PiecewiseCurve : time discontinuity between T_max_ and pol.min()");
        }
        curves_.push_back(cf);
        size_ = curves_.size();
        T_max_ = cf.max();
        time_curves_.push_back(T_max_);
    }

    ///  \brief Check if the curve is continuous of order given.
    ///  \param order : order of continuity we want to check.
    ///  \return True if the curve is continuous of order given.
    ///
    bool is_continuous(const std::size_t order)
    {
        bool isContinuous = true;
        std::size_t i=0;
        point_t value_end, value_start;
        while (isContinuous && i<(size_-1))
        {
            curve_t current = curves_.at(i);
            curve_t next = curves_.at(i+1);

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

            if ((value_end-value_start).norm() > MARGIN)
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
    std::size_t find_interval(const Numeric t) const
    {   
        // time before first control point time.
        if(t < time_curves_[0])
        {
            return 0;
        }
        // time is after last control point time
        if(t > time_curves_[size_-1])
        {
            return size_-1;
        }

        std::size_t left_id = 0;
        std::size_t right_id = size_-1;
        while(left_id <= right_id)
        {
            const std::size_t middle_id = left_id + (right_id - left_id)/2;
            if(time_curves_.at(middle_id) < t)
            {
                left_id = middle_id+1;
            }
            else if(time_curves_.at(middle_id) > t)
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
    t_curve_t curves_; // for curves 0/1/2 : [ curve0, curve1, curve2 ]
    t_time_t time_curves_; // for curves 0/1/2 : [ Tmin0, Tmax0,Tmax1,Tmax2 ]
    std::size_t size_; // Number of segments in piecewise curve = size of curves_
    Time T_min_, T_max_;
    static const double MARGIN;
};

template<typename Time, typename Numeric, std::size_t Dim, bool Safe, typename Point, typename T_Point, typename Curve>
const double piecewise_curve<Time, Numeric, Dim, Safe, Point, T_Point, Curve>::MARGIN(0.001);

} // end namespace


#endif // _CLASS_PIECEWISE_CURVE