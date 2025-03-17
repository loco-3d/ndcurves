/**
 * \file bezier_curve.h
 * \brief class allowing to create a Bezier curve of dimension 1 <= n <= 3.
 * \author Steve T.
 * \version 0.1
 * \date 06/17/2013
 */

#ifndef _CLASS_BERNSTEIN
#define _CLASS_BERNSTEIN

#include <math.h>

#include <stdexcept>
#include <vector>

#include "MathDefs.h"
#include "curve_abc.h"

namespace ndcurves {
/// \brief Computes a binomial coefficient  .
/// \param n : an unsigned integer.
/// \param k : an unsigned integer.
/// \return \f$\binom{n}{k}f$
///
inline unsigned int bin(const unsigned int n, const unsigned int k) {
  if (k > n)
    throw std::runtime_error("binomial coefficient higher than degree");
  if (k == 0) return 1;
  if (k > n / 2) return bin(n, n - k);
  return n * bin(n - 1, k - 1) / k;
}

/// \class Bernstein.
/// \brief Computes a Bernstein polynomial.
///
template <typename Numeric = double>
struct Bern {
  Bern() {}
  Bern(const unsigned int m, const unsigned int i)
      : m_minus_i(m - i), i_(i), bin_m_i_(bin(m, i)) {}

  virtual ~Bern() {}

  /// \brief Evaluation of Bernstein polynomial at value u.
  /// \param u : value between 0 and 1.
  /// \return Evaluation corresponding at value u.
  Numeric operator()(const Numeric u) const {
    if (!(u >= 0. && u <= 1.)) {
      throw std::invalid_argument("u needs to be betwen 0 and 1.");
    }
    return bin_m_i_ * (pow(u, i_)) * pow((1 - u), m_minus_i);
  }

  /// \brief Check if actual Bernstein polynomial and other are approximately
  /// equal. \param other : the other Bernstein polynomial to check. \return
  /// true if the two Bernstein polynomials are approximately equals.
  virtual bool operator==(const Bern& other) const {
    return ndcurves::isApprox<Numeric>(m_minus_i, other.m_minus_i) &&
           ndcurves::isApprox<Numeric>(i_, other.i_) &&
           ndcurves::isApprox<Numeric>(bin_m_i_, other.bin_m_i_);
  }

  /// \brief Check if actual Bernstein polynomial and other are different.
  /// \param other : the other Bernstein polynomial to check.
  /// \return true if the two Bernstein polynomials are different.
  virtual bool operator!=(const Bern& other) const { return !(*this == other); }

  /* Attributes */
  Numeric m_minus_i;
  Numeric i_;
  Numeric bin_m_i_;
  /* Attributes */

  // Serialization of the class
  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version) {
    if (version) {
      // Do something depending on version ?
    }
    ar& boost::serialization::make_nvp("m_minus_i", m_minus_i);
    ar& boost::serialization::make_nvp("i", i_);
    ar& boost::serialization::make_nvp("bin_m_i", bin_m_i_);
  }
};  // End struct Bern

/// \brief Computes all Bernstein polynomes for a certain degree.
///
template <typename Numeric>
std::vector<Bern<Numeric> > makeBernstein(const unsigned int n) {
  std::vector<Bern<Numeric> > res;
  for (unsigned int i = 0; i <= n; ++i) {
    res.push_back(Bern<Numeric>(n, i));
  }
  return res;
}
}  // namespace ndcurves

DEFINE_CLASS_TEMPLATE_VERSION(typename Numeric, ndcurves::Bern<Numeric>)

#endif  //_CLASS_BERNSTEIN
