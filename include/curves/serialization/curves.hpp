// Copyright (c) 2020, CNRS
// Authors: Pierre Fernbach <pierre.fernbach@laas.fr>

#ifndef CURVES_SERIALIZAION
#define CURVES_SERIALIZAION

/**
 * This file must be included by all classes containing
 * pointer of abstract curves as member and serializing them.
 * curves::serialization::register_types(ar) must be called before serializing them.
 */

#include <boost/shared_ptr.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include "registeration.hpp"
#include "curves/curve_abc.h"
#include "curves/so3_linear.h"
#include "curves/se3_curve.h"
#include "curves/polynomial.h"
#include "curves/bezier_curve.h"
#include "curves/piecewise_curve.h"
#include "curves/exact_cubic.h"
#include "curves/cubic_hermite_spline.h"


#endif  // ifndef CURVES_SERIALIZAION
