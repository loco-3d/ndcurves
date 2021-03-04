// Copyright (c) 2020, CNRS
// Authors: Pierre Fernbach <pierre.fernbach@laas.fr>

#ifndef CURVES_SERIALIZAION
#define CURVES_SERIALIZAION

/**
 * This file must be included by all classes containing
 * pointer of abstract curves as member and serializing them.
 * nd_curves::serialization::register_types(ar) must be called before serializing them.
 */

#include <boost/shared_ptr.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include "registeration.hpp"
#include "nd_curves/curve_abc.h"
#include "nd_curves/so3_linear.h"
#include "nd_curves/se3_curve.h"
#include "nd_curves/sinusoidal.h"
#include "nd_curves/polynomial.h"
#include "nd_curves/bezier_curve.h"
#include "nd_curves/constant_curve.h"
#include "nd_curves/piecewise_curve.h"
#include "nd_curves/exact_cubic.h"
#include "nd_curves/cubic_hermite_spline.h"


#endif  // ifndef CURVES_SERIALIZAION
