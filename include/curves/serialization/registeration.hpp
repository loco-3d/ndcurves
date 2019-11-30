/**
 * \file registeration.h
 * \brief registeration of class for serialization
 * \author Pierre Fernbach
 * \version 0.1
 * \date 27/11/19
 *
 * Boost::serialization need to be aware of all the derived class of the abstract curve_abc,
 * and thier template parameters.
 * Otherwise we cannot serialize a pointer to curve_abc.
 * New class should be added at the end and the order should not change
 * or it will break backward compatibility when deserializing objects.
 *
 */


#ifndef CURVES_REGISTERATION_H
#define CURVES_REGISTERATION_H
#include <Eigen/Dense>
#include <vector>

#include "curves/fwd.h"

namespace curves {
namespace serialization {

  template <class Archive>
  void register_types(Archive& ar){
    typedef double Scalar;
    typedef Eigen::Matrix<Scalar, -1, 1> pointX_t;
    typedef std::vector<pointX_t, Eigen::aligned_allocator<pointX_t> > t_point_t;

    //register derived class
    ar.template register_type<polynomial_t >();
    ar.template register_type<exact_cubic_t >();
    ar.template register_type<bezier_curve_t >();
    ar.template register_type<SO3Linear_t >();
    ar.template register_type<SE3Curve_t >();
    ar.template register_type<cubic_hermite_spline_t >();
    ar.template register_type<piecewise_curve_t >();
    ar.template register_type<piecewise_SE3_curve_t >();

  }

}
}

#endif // CURVES_REGISTERATION_H
