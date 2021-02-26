// Copyright (c) 2020, CNRS
// Authors: Pierre Fernbach <pierre.fernbach@laas.fr>

#ifndef CURVES_SERIALIZAION
#define CURVES_SERIALIZAION

/**
 * This file must be included by all classes containing
 * pointer of abstract curves as member and serializing them.
 * NDcurves::serialization::register_types(ar) must be called before serializing them.
 */

#include <boost/shared_ptr.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include "registeration.hpp"
#include "NDcurves/curve_abc.h"
#include "NDcurves/so3_linear.h"
#include "NDcurves/se3_curve.h"
#include "NDcurves/sinusoidal.h"
#include "NDcurves/polynomial.h"
#include "NDcurves/bezier_curve.h"
#include "NDcurves/constant_curve.h"
#include "NDcurves/piecewise_curve.h"
#include "NDcurves/exact_cubic.h"
#include "NDcurves/cubic_hermite_spline.h"


#endif  // ifndef CURVES_SERIALIZAION
