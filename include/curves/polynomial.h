/**
* \file polynomial.h
* \brief Definition of a cubic spline.
* \author Steve T.
* \version 0.1
* \date 06/17/2013
*
* This file contains definitions for the polynomial struct.
* It allows the creation and evaluation of natural
* smooth splines of arbitrary dimension and order
*/


#ifndef _STRUCT_POLYNOMIAL
#define _STRUCT_POLYNOMIAL

#include "MathDefs.h"

#include "curve_abc.h"

#include <iostream>
#include <algorithm>
#include <functional>
#include <stdexcept>

#include "serialization/archive.hpp"
#include "serialization/eigen-matrix.hpp"

namespace curves
{
/// \class polynomial.
/// \brief Represents a polynomial of an arbitrary order defined on the interval
/// \f$[t_{min}, t_{max}]\f$. It follows the equation :<br>
/// \f$ x(t) = a + b(t - t_{min}) + ... + d(t - t_{min})^N \f$<br> 
/// where N is the order and \f$ t \in [t_{min}, t_{max}] \f$.
///
template<typename Time= double, typename Numeric=Time, std::size_t Dim=3, bool Safe=false,
         typename Point= Eigen::Matrix<Numeric, Eigen::Dynamic, 1>, typename T_Point =std::vector<Point,Eigen::aligned_allocator<Point> > >
struct polynomial : public curve_abc<Time, Numeric, Safe, Point>,
                    public serialization::Serializable< polynomial<Time, Numeric, Dim, Safe, Point, T_Point> >
{
    typedef Point 	point_t;
    typedef T_Point t_point_t;
    typedef Time 	time_t;
    typedef Numeric	num_t;
    typedef curve_abc<Time, Numeric, Safe, Point> curve_abc_t;
    typedef Eigen::Matrix<double, Dim, Eigen::Dynamic> coeff_t;
    typedef Eigen::Ref<coeff_t> coeff_t_ref;

/* Constructors - destructors */
    public:

    polynomial() {}

    /// \brief Constructor.
    /// \param coefficients : a reference to an Eigen matrix where each column is a coefficient,
    /// from the zero order coefficient, up to the highest order. Spline order is given
    /// by the number of the columns -1.
    /// \param min  : LOWER bound on interval definition of the curve.
    /// \param max  : UPPER bound on interval definition of the curve.
    polynomial(const coeff_t& coefficients, const time_t min, const time_t max)
        : curve_abc_t(),
          coefficients_(coefficients), 
          dim_(Dim), degree_(coefficients_.cols()-1), 
          T_min_(min), T_max_(max)
    {
        safe_check();
    }

    /// \brief Constructor
    /// \param coefficients : a container containing all coefficients of the spline, starting
    ///  with the zero order coefficient, up to the highest order. Spline order is given
    ///  by the size of the coefficients.
    /// \param min  : LOWER bound on interval definition of the spline.
    /// \param max  : UPPER bound on interval definition of the spline.
    polynomial(const T_Point& coefficients, const time_t min, const time_t max)
        : curve_abc_t(),
          coefficients_(init_coeffs(coefficients.begin(), coefficients.end())),
          dim_(Dim), degree_(coefficients_.cols()-1), 
          T_min_(min), T_max_(max)
    {
        safe_check();
    }

    /// \brief Constructor.
    /// \param zeroOrderCoefficient : an iterator pointing to the first element of a structure containing the coefficients
    ///  it corresponds to the zero degree coefficient.
    /// \param out   : an iterator pointing to the last element of a structure ofcoefficients.
    /// \param min   : LOWER bound on interval definition of the spline.
    /// \param max   : UPPER bound on interval definition of the spline.
    template<typename In>
    polynomial(In zeroOrderCoefficient, In out, const time_t min, const time_t max)
        : coefficients_(init_coeffs(zeroOrderCoefficient, out)),
          dim_(Dim), degree_(coefficients_.cols()-1), 
          T_min_(min), T_max_(max)
    {
        safe_check();
    }

    /// \brief Destructor
    ~polynomial()
    {
        // NOTHING
    }


    polynomial(const polynomial& other)
        : coefficients_(other.coefficients_),
          dim_(other.dim_), degree_(other.degree_), T_min_(other.T_min_), T_max_(other.T_max_)
    {}


    //polynomial& operator=(const polynomial& other);

    private:
    void safe_check()
    {
        if(Safe)
        {
            if(T_min_ > T_max_)
            {
                std::invalid_argument("Tmin should be inferior to Tmax");
            }
            if(coefficients_.size() != int(degree_+1))
            {
                std::runtime_error("Spline order and coefficients do not match");
            }
        }
    }

/* Constructors - destructors */

/*Operations*/
    public:
    /*///  \brief Evaluation of the cubic spline at time t.
    ///  \param t : the time when to evaluate the spine
    ///  \param return \f$x(t)\f$, point corresponding on curve at time t.
    virtual point_t operator()(const time_t t) const
    {
        if((t < T_min_ || t > T_max_) && Safe){ throw std::out_of_range("TODO");}
        time_t const dt (t-T_min_);
        time_t cdt(1);
        point_t currentPoint_ = point_t::Zero();
        for(int i = 0; i < degree_+1; ++i, cdt*=dt)
            currentPoint_ += cdt *coefficients_.col(i);
        return currentPoint_;
    }*/


    ///  \brief Evaluation of the cubic spline at time t using horner's scheme.
    ///  \param t : time when to evaluate the spline.
    ///  \return \f$x(t)\f$ point corresponding on spline at time t.
    virtual point_t operator()(const time_t t) const
    {
        if((t < T_min_ || t > T_max_) && Safe)
        { 
            throw std::invalid_argument("error in polynomial : time t to evaluate should be in range [Tmin, Tmax] of the curve");
        }
        time_t const dt (t-T_min_);
        point_t h = coefficients_.col(degree_);
        for(int i=(int)(degree_-1); i>=0; i--)
        {
            h = dt*h + coefficients_.col(i);
        }
        return h;
    }


    ///  \brief Evaluation of the derivative of order N of spline at time t.
    ///  \param t : the time when to evaluate the spline.
    ///  \param order : order of derivative.
    ///  \return \f$\frac{d^Nx(t)}{dt^N}\f$ point corresponding on derivative spline at time t.
    virtual point_t derivate(const time_t t, const std::size_t order) const
    {
        if((t < T_min_ || t > T_max_) && Safe)
        { 
            throw std::invalid_argument("error in polynomial : time t to evaluate derivative should be in range [Tmin, Tmax] of the curve");
        }
        time_t const dt (t-T_min_);
        time_t cdt(1);
        point_t currentPoint_ = point_t::Zero(dim_);
        for(int i = (int)(order); i < (int)(degree_+1); ++i, cdt*=dt) 
        {
            currentPoint_ += cdt *coefficients_.col(i) * fact(i, order);
        }
        return currentPoint_;
    }

    private:
    num_t fact(const std::size_t n, const std::size_t order) const
    {
        num_t res(1);
        for(std::size_t i = 0; i < order; ++i)
        {
            res *= (num_t)(n-i);
        }
        return res;
    }

/*Operations*/

/*Helpers*/
    public:
    /// \brief Get the minimum time for which the curve is defined
    /// \return \f$t_{min}\f$ lower bound of time range.
    num_t virtual min() const {return T_min_;}
    /// \brief Get the maximum time for which the curve is defined.
    /// \return \f$t_{max}\f$ upper bound of time range.
    num_t virtual max() const {return T_max_;}
/*Helpers*/

/*Attributes*/
    public:
    coeff_t coefficients_; //const
    std::size_t dim_; //const
    std::size_t degree_; //const

    private:
    time_t T_min_, T_max_;
/*Attributes*/

    private:

    template<typename In>
    coeff_t init_coeffs(In zeroOrderCoefficient, In highestOrderCoefficient)
    {
        std::size_t size = std::distance(zeroOrderCoefficient, highestOrderCoefficient);
        coeff_t res = coeff_t(Dim, size); int i = 0;
        for(In cit = zeroOrderCoefficient; cit != highestOrderCoefficient; ++cit, ++i)
        {
            res.col(i) = *cit;
        }
        return res;
    }

    public:
    // Serialization of the class
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive& ar, const unsigned int version){
        ar & boost::serialization::make_nvp("coefficients", coefficients_);
        ar & boost::serialization::make_nvp("dim", dim_);
        ar & boost::serialization::make_nvp("degree", degree_);
        ar & boost::serialization::make_nvp("T_min", T_min_);
        ar & boost::serialization::make_nvp("T_max", T_max_);
    }

}; //class polynomial

} // namespace curves
#endif //_STRUCT_POLYNOMIAL

