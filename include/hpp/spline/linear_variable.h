/**
* \file linear_variable.h
* \brief storage for variable points of the form p_i = a_i x + b_i
* \author Steve T.
* \version 0.1
* \date 07/02/2019
*/


#ifndef _CLASS_LINEAR_VARIABLE
#define _CLASS_LINEAR_VARIABLE

#include "curve_abc.h"

#include "MathDefs.h"

#include <math.h>
#include <vector>
#include <Eigen/Core>
#include <stdexcept>

namespace spline
{
template <int Dim, typename Numeric=double>
struct linear_variable{
    typedef Numeric num_t;
    typedef Eigen::Matrix<num_t, Dim, Dim> matrix_t;
    typedef Eigen::Matrix<num_t, Dim, 1> point_t;

    typedef linear_variable<Dim, Numeric> linear_variable_t;

    matrix_t A_;
    point_t b_;

    linear_variable(): A_(matrix_t::Identity()), b_(point_t::Zero()){}
    linear_variable(const matrix_t& A, const point_t& b):A_(A),b_(b) {}
    linear_variable(const point_t& b):A_(matrix_t::Zero()),b_(b) {} // constant


    linear_variable& operator+=(const linear_variable& w1)
    {
        this->A_ += w1.A_;
        this->b_ += w1.b_;
        return *this;
    }
    linear_variable& operator-=(const linear_variable& w1)
    {
        this->A_ -= w1.A_;
        this->b_ -= w1.b_;
        return *this;
    }

    static linear_variable_t Zero(size_t dim=0){
        linear_variable_t w;
        w.A_  = matrix_t::Zero();
        w.b_  = point_t::Zero();
        return w;
    }
};


template<typename Var>
struct variables{
    typedef Var var_t;
    typedef variables<Var> variables_t;

    typedef std::vector<var_t> T_var_t;
    typedef typename T_var_t::iterator  IT_var_t;
    typedef typename T_var_t::const_iterator CIT_var_t;

    T_var_t variables_;

    variables() {}

    variables& operator+=(const variables& w1)
    {
        if(variables_.size() == 0)
            variables_ = w1.variables_;
        else if (w1.variables_.size() !=0)
        {
            assert(variables_.size() == w1.variables_.size());
            CIT_var_t cit = w1.variables_.begin();
            for(IT_var_t it = variables_.begin(); it != variables_.end(); ++it, ++cit)
                (*it)+=(*cit);
        }
        return *this;
    }

    variables& operator-=(const variables& w1)
    {
        if(variables_.size() == 0)
            variables_ = w1.variables_;
        else if (w1.variables_.size() !=0)
        {
            assert(variables_.size() == w1.variables_.size());
            CIT_var_t cit = w1.variables_.begin();
            for(IT_var_t it = variables_.begin(); it != variables_.end(); ++it, ++cit)
                (*it)-=(*cit);
        }
        return *this;
    }

    static variables_t Zero(size_t /*dim*/){
        variables_t w;
        return w;
    }
};

template <int D, typename N>
inline linear_variable<D,N> operator+(const linear_variable<D,N>& w1, const linear_variable<D,N>& w2)
{
    return linear_variable<D,N>(w1.A_ + w2.A_, w1.b_ + w2.b_);
}

template <int D, typename N>
linear_variable<D,N> operator-(const linear_variable<D,N>& w1, const linear_variable<D,N>& w2)
{
    return linear_variable<D,N>(w1.A_ - w2.A_, w1.b_ - w2.b_);
}

template <int D, typename N>
linear_variable<D,N> operator*(const double k, const linear_variable<D,N>& w){
    return linear_variable<D,N>(k*w.A_,k*w.b_);
}

template <int D, typename N>
linear_variable<D,N> operator*(const linear_variable<D,N>& w,const double k){
    return linear_variable<D,N>(k*w.A_,k*w.b_);
}

template <int D, typename N>
linear_variable<D,N> operator/(const linear_variable<D,N>& w,const double k){
    return linear_variable<D,N>(w.A_/k,w.b_/k);
}

template<typename V>
variables<V> operator+(const variables<V>& w1, const variables<V>& w2)
{
    if(w2.variables_.size() == 0)
        return w1;
    if(w1.variables_.size() == 0)
        return w2;
    variables<V> res;
    assert(w2.variables_.size() == w1.variables_.size());
    typename variables<V>::CIT_var_t cit = w1.variables_.begin();
    for(typename variables<V>::CIT_var_t cit2 = w2.variables_.begin(); cit2 != w2.variables_.end(); ++cit, ++cit2)
        res.variables_.push_back((*cit)+(*cit2));
    return res;
}

template<typename V>
variables<V> operator-(const variables<V>& w1, const variables<V>& w2)
{
    if(w2.variables_.size() == 0)
        return w1;
    if(w1.variables_.size() == 0)
        return w2;
    variables<V> res;
    assert(w2.variables_.size() == w1.variables_.size());
    typename variables<V>::CIT_var_t cit = w1.variables_.begin();
    for(typename variables<V>::CIT_var_t cit2 = w2.variables_.begin(); cit2 != w2.variables_.end(); ++cit, ++cit2)
        res.variables_.push_back((*cit)-(*cit2));
    return res;
}

template<typename V>
variables<V> operator*(const double k, const variables<V>& w)
{
    if(w.variables_.size() == 0)
        return w;
    variables<V> res;
    for(typename variables<V>::CIT_var_t cit = w.variables_.begin(); cit != w.variables_.end(); ++cit)
        res.variables_.push_back(k*(*cit));
    return res;
}

template<typename V>
variables<V> operator*(const variables<V>& w,const double k)
{
    if(w.variables_.size() == 0)
        return w;
    variables<V> res;
    for(typename variables<V>::CIT_var_t cit = w.variables_.begin(); cit != w.variables_.end(); ++cit)
        res.variables_.push_back((*cit)*k);
    return res;
}

template<typename V>
variables<V> operator/(const variables<V>& w,const double k)
{
    if(w.variables_.size() == 0)
        return w;
    variables<V> res;
    for(typename variables<V>::CIT_var_t cit = w.variables_.begin(); cit != w.variables_.end(); ++cit)
        res.variables_.push_back((*cit)/k);
    return res;
}

} // namespace spline
#endif //_CLASS_LINEAR_VARIABLE

