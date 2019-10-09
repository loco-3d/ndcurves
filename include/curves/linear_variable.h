/**
* \file linear_variable.h
* \brief storage for variable points of the form p_i = B_i x + c_i
* \author Steve T.
* \version 0.1
* \date 07/02/2019
*/


#ifndef _CLASS_LINEAR_VARIABLE
#define _CLASS_LINEAR_VARIABLE

#include "curve_abc.h"
#include "bezier_curve.h"

#include "MathDefs.h"

#include <math.h>
#include <vector>
#include <Eigen/Core>
#include <stdexcept>

namespace curves
{
template <typename Numeric=double, bool Safe=true>
struct linear_variable
{
    typedef Eigen::Matrix<Numeric, Eigen::Dynamic, Eigen::Dynamic> matrix_dim_t;
    typedef Eigen::Matrix<Numeric, Eigen::Dynamic, Eigen::Dynamic> matrix_dim_x_t;
    typedef Eigen::Matrix<Numeric, Eigen::Dynamic, 1> point_dim_t;
    typedef Eigen::Matrix<Numeric, Eigen::Dynamic, 1> vector_x_t;
    typedef Eigen::Matrix<Numeric, Eigen::Dynamic, Eigen::Dynamic> matrix_x_t;
    typedef linear_variable<Numeric> linear_variable_t;

    linear_variable(): B_(matrix_dim_t::Identity(0,0)), c_(point_dim_t::Zero(0)), zero(true){} //variable
    linear_variable(const point_dim_t& c):B_(matrix_dim_t::Zero(c.size(),c.size())),c_(c), zero(false) {} // constant
    linear_variable(const matrix_dim_x_t& B, const point_dim_t& c):B_(B),c_(c), zero(false) {} //mixed

    // linear evaluation
    point_dim_t operator()(const Eigen::Ref<const vector_x_t>& val) const
    {
        if(isZero())
            return c();
        if(Safe && B().cols() != val.rows() )
            throw std::length_error("Cannot evaluate linear variable, variable value does not have the correct dimension");
        return B() * val + c();
    }

    linear_variable_t& operator+=(const linear_variable_t& w1)
    {
        if (w1.isZero())
            return *this;
        if(isZero())
        {
            this->B_ = w1.B_;
            zero = w1.isZero();
        }
        else
        {
            this->B_ += w1.B_;
        }
        this->c_ += w1.c_;
        return *this;
    }
    linear_variable_t& operator-=(const linear_variable_t& w1)
    {
        if (w1.isZero())
            return *this;
        if(isZero())
        {
            this->B_ = -w1.B_;
            zero = w1.isZero();
        }
        else
        {
            this->B_ -= w1.B_;
        }
        this->c_ -= w1.c_;
        return *this;
    }
    linear_variable_t& operator/=(const double d)
    {
        B_ /= d;
        c_ /= d;
        return *this;
    }
    linear_variable_t& operator*=(const double d)
    {
        B_ *= d;
        c_ *= d;
        return *this;
    }

    static linear_variable_t Zero(size_t dim=0){
        return linear_variable_t(matrix_x_t::Identity(dim,dim), vector_x_t::Zero(dim));
    }

    std::size_t size() const {return zero ? 0 : std::max(B_.cols(), c_.size()) ;}

    const matrix_dim_x_t& B() const {return B_;}
    const point_dim_t& c () const {return c_;}
    bool isZero () const {return zero;}

private:
    matrix_dim_x_t B_;
    point_dim_t c_;
    bool zero;
};

template <typename N, bool S>
inline linear_variable<N,S> operator+(const linear_variable<N,S>& w1, const linear_variable<N,S>& w2)
{
    linear_variable<N,S> res(w1.B(), w1.c());
    return res+=w2;
}

template <typename N, bool S>
linear_variable<N,S> operator-(const linear_variable<N,S>& w1, const linear_variable<N,S>& w2)
{
    linear_variable<N,S> res(w1.B(), w1.c());
    return res-=w2;
}

template <typename N, bool S>
linear_variable<N,S> operator*(const double k, const linear_variable<N,S>& w){
    linear_variable<N,S> res(w.B(), w.c());
    return res*=k;
}

template <typename N, bool S>
linear_variable<N,S> operator*(const linear_variable<N,S>& w,const double k)
{
    linear_variable<N,S> res(w.B(), w.c());
    return res*=k;
}

template <typename N, bool S>
linear_variable<N,S> operator/(const linear_variable<N,S>& w,const double k)
{
    linear_variable<N,S> res(w.B(), w.c());
    return res/=k;
}

template <typename BezierFixed, typename BezierLinear, typename X >
BezierFixed evaluateLinear(const BezierLinear& bIn, const X x)
{
    typename BezierFixed::t_point_t fixed_wps;
    for (typename BezierLinear::cit_point_t cit = bIn.waypoints().begin(); cit != bIn.waypoints().end(); ++cit)
        fixed_wps.push_back(cit->operator()(x));
    return BezierFixed(fixed_wps.begin(),fixed_wps.end(), bIn.T_min_,bIn.T_max_);
}

} // namespace curves
#endif //_CLASS_LINEAR_VARIABLE

