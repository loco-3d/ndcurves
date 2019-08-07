/**
* \file bezier_curve.h
* \brief class allowing to create a Bezier curve of dimension 1 <= n <= 3.
* \author Steve T.
* \version 0.1
* \date 06/17/2013
*/


#ifndef _CLASS_BERNSTEIN
#define _CLASS_BERNSTEIN

#include "curve_abc.h"

#include "MathDefs.h"

#include "serialization/archive.hpp"

#include <math.h>
#include <vector>
#include <stdexcept>

namespace curves
{
/// \brief Computes a binomial coefficient.
/// \param n : an unsigned integer.
/// \param k : an unsigned integer.
/// \return \f$\binom{n}{k}f$
///
inline unsigned int bin(const unsigned  int n, const unsigned  int k)
{
    if(k >  n)
        throw std::runtime_error("binomial coefficient higher than degree");
    if(k == 0)
        return 1;
    if(k > n/2)
        return bin(n,n-k);
    return n * bin(n-1,k-1) / k;
}

/// \class Bernstein.
/// \brief Computes a Bernstein polynome.
///
template <typename Numeric = double>
struct Bern : public serialization::Serializable< Bern<Numeric> > {

Bern(){}

Bern(const unsigned int m, const unsigned int i)
    :m_minus_i(m - i)
    ,i_(i)
    ,bin_m_i_(bin(m,i)) {}

~Bern(){}

Numeric operator()(const Numeric u) const
{
    assert(u >= 0. && u <= 1.);
    return bin_m_i_*(pow(u, i_)) *pow((1-u),m_minus_i);
}

Numeric m_minus_i;
Numeric i_;
Numeric bin_m_i_;

// Serialization of the class
friend class boost::serialization::access;

template<class Archive>
void serialize(Archive& ar, const unsigned int version){
    ar & boost::serialization::make_nvp("m_minus_i", m_minus_i);
    ar & boost::serialization::make_nvp("i", i_);
    ar & boost::serialization::make_nvp("bin_m_i", bin_m_i_);
}
};

/// \brief Computes all Bernstein polynomes for a certain degree.
///
template <typename Numeric>
std::vector<Bern<Numeric> > makeBernstein(const unsigned int n)
{
    std::vector<Bern<Numeric> > res;
    for(unsigned int i = 0; i<= n; ++i)
        res.push_back(Bern<Numeric>(n, i));
    return res;
}
} // namespace curves
#endif //_CLASS_BERNSTEIN

