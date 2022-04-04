#include <Eigen/Dense>
#include <eigenpy/eigenpy.hpp>
#include <eigenpy/geometry.hpp>
#include <eigenpy/memory.hpp>
#include <vector>

#include "ndcurves/bernstein.h"
#include "ndcurves/bezier_curve.h"
#include "ndcurves/constant_curve.h"
#include "ndcurves/cubic_hermite_spline.h"
#include "ndcurves/curve_constraint.h"
#include "ndcurves/curve_conversion.h"
#include "ndcurves/exact_cubic.h"
#include "ndcurves/fwd.h"
#include "ndcurves/linear_variable.h"
#include "ndcurves/piecewise_curve.h"
#include "ndcurves/polynomial.h"
#include "ndcurves/python/python_definitions.h"
#include "ndcurves/se3_curve.h"
#include "ndcurves/sinusoidal.h"
#include "ndcurves/so3_linear.h"

#ifndef _VARIABLES_PYTHON_BINDINGS
#define _VARIABLES_PYTHON_BINDINGS

namespace ndcurves {
static const int dim = 3;
typedef linear_variable<real, true> linear_variable_t;
typedef quadratic_variable<real> quadratic_variable_t;
typedef bezier_curve<real, real, true, linear_variable_t>
    bezier_linear_variable_t;

/*linear variable control points*/
bezier_linear_variable_t* wrapBezierLinearConstructor(
    const point_list3_t& matrices, const point_list3_t& vectors);
bezier_linear_variable_t* wrapBezierLinearConstructorBounds(
    const point_list3_t& matrices, const point_list3_t& vectors,
    const real T_min, const real T_max);

typedef Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic> matrix_x_t;

struct matrix_pair {
  matrix_pair() {}
  matrix_pair(const Eigen::Ref<const matrix_x_t> A,
              const Eigen::Ref<const matrix_x_t> b)
      : A_(A), b_(b) {}
  matrix_x_t A_;
  matrix_x_t b_;
  matrix_x_t A() { return A_; }
  matrix_x_t b() { return b_; }
};

Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic> cost_t_quad(
    const quadratic_variable_t& p);
Eigen::Matrix<real, Eigen::Dynamic, 1> cost_t_linear(
    const quadratic_variable_t& p);
real cost_t_constant(const quadratic_variable_t& p);

linear_variable_t* wayPointsToLists(const bezier_linear_variable_t& self);

struct LinearBezierVector {
  std::vector<bezier_linear_variable_t> beziers_;
  std::size_t size() { return beziers_.size(); }
  bezier_linear_variable_t* at(std::size_t i) {
    assert(i < size());
    return new bezier_linear_variable_t(beziers_[i]);
  }
};

/*** TEMPLATE SPECIALIZATION FOR PYTHON ****/
}  // namespace ndcurves
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(ndcurves::bernstein_t)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(ndcurves::curve_constraints_t)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(ndcurves::matrix_x_t)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(ndcurves::pointX_t)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(ndcurves::linear_variable_t)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(
    ndcurves::bezier_linear_variable_t)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(ndcurves::matrix_pair)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(ndcurves::polynomial_t)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(ndcurves::exact_cubic_t)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(ndcurves::bezier_t)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(ndcurves::constant_t)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(ndcurves::cubic_hermite_spline_t)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(ndcurves::piecewise_t)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(ndcurves::polynomial3_t)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(ndcurves::exact_cubic3_t)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(ndcurves::bezier3_t)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(ndcurves::constant3_t)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(
    ndcurves::cubic_hermite_spline3_t)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(ndcurves::piecewise3_t)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(ndcurves::SO3Linear_t)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(ndcurves::SE3Curve_t)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(ndcurves::sinusoidal_t)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(ndcurves::piecewise_SE3_t)

#endif  //_VARIABLES_PYTHON_BINDINGS
