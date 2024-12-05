// Copyright (c) 2020, CNRS
// Authors: Pierre Fernbach <pierre.fernbach@laas.fr>

#ifndef CURVES_SERIALIZAION
#define CURVES_SERIALIZAION

/**
 * This file must be included by all classes containing
 * pointer of abstract curves as member and serializing them.
 * ndcurves::serialization::register_types(ar) must be called before serializing
 * them.
 */

#include <boost/serialization/shared_ptr.hpp>

#include "ndcurves/bezier_curve.h"
#include "ndcurves/constant_curve.h"
#include "ndcurves/cubic_hermite_spline.h"
#include "ndcurves/curve_abc.h"
#include "ndcurves/exact_cubic.h"
#include "ndcurves/piecewise_curve.h"
#include "ndcurves/polynomial.h"
#include "ndcurves/se3_curve.h"
#include "ndcurves/sinusoidal.h"
#include "ndcurves/so3_linear.h"
#include "ndcurves/so3_smooth.h"
#include "registeration.hpp"

#endif  // ifndef CURVES_SERIALIZAION
