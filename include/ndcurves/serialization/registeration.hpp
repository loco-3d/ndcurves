/**
 * \file registeration.h
 * \brief registeration of class for serialization
 * \author Pierre Fernbach
 * \version 0.1
 * \date 27/11/19
 *
 * Boost::serialization need to be aware of all the derived class of the
 * abstract curve_abc, and thier template parameters. Otherwise we cannot
 * serialize a pointer to curve_abc. New class should be added at the end and
 * the order should not change or it will break backward compatibility when
 * deserializing objects.
 *
 */

#ifndef CURVES_REGISTERATION_H
#define CURVES_REGISTERATION_H
#include <Eigen/Dense>
#include <vector>

#include "ndcurves/fwd.h"

/**
 *This file define a method register_types that
 * register all the curves class of this package for a boost::Archive
 * This is used to serialize pointer of the abstract class curve_abc
 */

namespace ndcurves {
namespace serialization {

template <class Archive>
void register_types(Archive& ar, const unsigned int version) {
  // register derived class
  ar.template register_type<polynomial_t>();
  ar.template register_type<exact_cubic_t>();
  ar.template register_type<bezier_t>();
  ar.template register_type<cubic_hermite_spline_t>();
  ar.template register_type<piecewise_t>();

  ar.template register_type<polynomial3_t>();
  ar.template register_type<exact_cubic3_t>();
  ar.template register_type<bezier3_t>();
  ar.template register_type<cubic_hermite_spline3_t>();
  ar.template register_type<piecewise3_t>();

  ar.template register_type<SO3Linear_t>();
  ar.template register_type<SE3Curve_t>();
  ar.template register_type<piecewise_SE3_t>();

  if (version >= 1) {
    ar.template register_type<constant3_t>();
    ar.template register_type<sinusoidal_t>();
    ar.template register_type<constant_t>();
    ar.template register_type<polynomial1_t>();
#ifdef CURVES_WITH_PINOCCHIO_SUPPORT
    ar.template register_type<SO3Smooth_t>();
#endif
  }
}

}  // namespace serialization
}  // namespace ndcurves

#endif  // CURVES_REGISTERATION_H
