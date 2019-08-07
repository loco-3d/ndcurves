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

#include "serialization/archive.hpp"
#include "serialization/eigen-matrix.hpp"

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
struct piecewise_curve : public curve_abc<Time, Numeric, Safe, Point>,
                         public serialization::Serializable< piecewise_curve<Time, Numeric, Dim, Safe, Point, T_Point, Curve> >
{
    typedef Point   point_t;
    typedef T_Point t_point_t;
    typedef Time    time_t;
    typedef Numeric num_t;
    typedef Curve curve_t;
    typedef typename std::vector < curve_t > t_curve_t;
    typedef typename std::vector< Time > t_time_t;

    public:

    piecewise_curve()
        : size_(0)
    {}

    /// \brief Constructor.
    /// Initialize a piecewise curve by giving the first curve.
    /// \param pol   : a polynomial curve.
    ///
    piecewise_curve(const curve_t& cf)
    {
        size_ = 0;
        add_curve(cf);
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
    ///         is equal to \f$T_{max}\f$ of the actual piecewise curve. The curve added should be of type Curve as defined
    ///         in the template.
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
        if (size_ == 1)
        {
            // First curve added
            time_curves_.push_back(cf.min());
            T_min_ = cf.min();
        }
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

    template<typename Bezier>
    piecewise_curve<Time, Numeric, Dim, Safe, Point, T_Point, Bezier> convert_piecewise_curve_to_bezier()
    {
        typedef piecewise_curve<Time, Numeric, Dim, Safe, Point, T_Point, Bezier> piecewise_curve_out_t;
        // Get first curve (segment)
        curve_t first_curve = curves_.at(0);
        Bezier first_curve_output = bezier_from_curve<Bezier, curve_t>(first_curve);
        // Create piecewise curve
        piecewise_curve_out_t pc_res(first_curve_output);
        // Convert and add all other curves (segments)
        for (std::size_t i=1; i<size_; i++)
        {
            pc_res.add_curve(bezier_from_curve<Bezier, curve_t>(curves_.at(i)));
        }
        return pc_res;
    }

    template<typename Hermite>
    piecewise_curve<Time, Numeric, Dim, Safe, Point, T_Point, Hermite> convert_piecewise_curve_to_cubic_hermite()
    {
        typedef piecewise_curve<Time, Numeric, Dim, Safe, Point, T_Point, Hermite> piecewise_curve_out_t;
        // Get first curve (segment)
        curve_t first_curve = curves_.at(0);
        Hermite first_curve_output = hermite_from_curve<Hermite, curve_t>(first_curve);
        // Create piecewise curve
        piecewise_curve_out_t pc_res(first_curve_output);
        // Convert and add all other curves (segments)
        for (std::size_t i=1; i<size_; i++)
        {
            pc_res.add_curve(hermite_from_curve<Hermite, curve_t>(curves_.at(i)));
        }
        return pc_res;
    }

    template<typename Polynomial>
    piecewise_curve<Time, Numeric, Dim, Safe, Point, T_Point, Polynomial> convert_piecewise_curve_to_polynomial()
    {
        typedef piecewise_curve<Time, Numeric, Dim, Safe, Point, T_Point, Polynomial> piecewise_curve_out_t;
        // Get first curve (segment)
        curve_t first_curve = curves_.at(0);
        Polynomial first_curve_output = polynomial_from_curve<Polynomial, curve_t>(first_curve);
        // Create piecewise curve
        piecewise_curve_out_t pc_res(first_curve_output);
        // Convert and add all other curves (segments)
        for (std::size_t i=1; i<size_; i++)
        {
            pc_res.add_curve(polynomial_from_curve<Polynomial, curve_t>(curves_.at(i)));
        }
        return pc_res;
    }

    template<typename Polynomial>
    static piecewise_curve<Time, Numeric, Dim, Safe, Point, T_Point, Polynomial> 
    convert_discrete_points_to_polynomial(T_Point points, Time T_min, Time T_max)
    {
        if(Safe &! (points.size()>1))
        {
            //std::cout<<"[Min,Max]=["<<T_min_<<","<<T_max_<<"]"<<" t="<<t<<std::endl;
            throw std::invalid_argument("piecewise_curve -> convert_discrete_points_to_polynomial, Error, less than 2 discrete points");
        }
        typedef piecewise_curve<Time, Numeric, Dim, Safe, Point, T_Point, Polynomial> piecewise_curve_out_t;
        Time discretization_step = (T_max-T_min)/Time(points.size()-1);
        Time time_actual = T_min;
        // Initialization at first points
        point_t actual_point = points[0];
        point_t next_point = points[1];
        point_t coeff_order_zero(actual_point);
        point_t coeff_order_one((next_point-actual_point)/discretization_step);
        t_point_t coeffs;
        coeffs.push_back(coeff_order_zero);
        coeffs.push_back(coeff_order_one);
        Polynomial pol(coeffs,time_actual,time_actual+discretization_step);
        piecewise_curve_out_t ppc(pol);
        time_actual += discretization_step;
        // Other points
        for (std::size_t i=1; i<points.size()-2; i++)
        {
            coeffs.clear();
            actual_point = points[i];
            next_point = points[i+1];
            coeff_order_zero = actual_point;
            coeff_order_one = (next_point-actual_point)/discretization_step;
            coeffs.push_back(coeff_order_zero);
            coeffs.push_back(coeff_order_one);
            ppc.add_curve(Polynomial(coeffs,time_actual,time_actual+discretization_step));
            time_actual += discretization_step;
        }
        // Last points
        coeffs.clear();
        actual_point = points[points.size()-2];
        next_point = points[points.size()-1];
        coeff_order_zero = actual_point;
        coeff_order_one = (next_point-actual_point)/discretization_step;
        coeffs.push_back(coeff_order_zero);
        coeffs.push_back(coeff_order_one);
        ppc.add_curve(Polynomial(coeffs,time_actual,T_max));
        return ppc;
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

    public:
    // Serialization of the class
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive& ar, const unsigned int version){
        ar & curves_;
        ar & time_curves_;
        ar & size_;
        ar & T_min_ & T_max_;
    }

};

template<typename Time, typename Numeric, std::size_t Dim, bool Safe, typename Point, typename T_Point, typename Curve>
const double piecewise_curve<Time, Numeric, Dim, Safe, Point, T_Point, Curve>::MARGIN(0.001);

} // end namespace


#endif // _CLASS_PIECEWISE_CURVE