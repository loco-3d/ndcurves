#include "curves/fwd.h"
#include "curves/linear_variable.h"
#include "curves/bezier_curve.h"
#include "curves/polynomial.h"
#include "curves/exact_cubic.h"
#include "curves/curve_constraint.h"
#include "curves/curve_conversion.h"
#include "curves/bernstein.h"
#include "curves/cubic_hermite_spline.h"
#include "curves/piecewise_curve.h"
#include "curves/so3_linear.h"
#include "curves/se3_curve.h"
#include "python_definitions.h"
#include <eigenpy/memory.hpp>
#include <eigenpy/eigenpy.hpp>
#include <eigenpy/geometry.hpp>
#include <Eigen/Dense>

#include <vector>

#ifndef _VARIABLES_PYTHON_BINDINGS
#define _VARIABLES_PYTHON_BINDINGS

namespace curves {
static const int dim = 3;
typedef linear_variable<real, true> linear_variable_t;
typedef quadratic_variable<real> quadratic_variable_t;
typedef bezier_curve<real, real, true, linear_variable_t> bezier_linear_variable_t;

/*linear variable control points*/
bezier_linear_variable_t* wrapBezierLinearConstructor(const point_list3_t& matrices, const point_list3_t& vectors);
bezier_linear_variable_t* wrapBezierLinearConstructorBounds(const point_list3_t& matrices, const point_list3_t& vectors,
                                                            const real T_min, const real T_max);

typedef Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic> matrix_x_t;

struct matrix_pair {
  matrix_pair() {}
  matrix_pair(const Eigen::Ref<const matrix_x_t> A, const Eigen::Ref<const matrix_x_t> b) : A_(A), b_(b) {}
  matrix_x_t A_;
  matrix_x_t b_;
  matrix_x_t A() { return A_; }
  matrix_x_t b() { return b_; }
};

Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic> cost_t_quad(const quadratic_variable_t& p);
Eigen::Matrix<real, Eigen::Dynamic, 1> cost_t_linear(const quadratic_variable_t& p);
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
}  // namespace curves
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(curves::bernstein_t)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(curves::curve_constraints_t)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(curves::matrix_x_t)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(curves::pointX_t)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(curves::linear_variable_t)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(curves::bezier_linear_variable_t)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(curves::matrix_pair)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(curves::polynomial_t)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(curves::exact_cubic_t)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(curves::bezier_t)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(curves::cubic_hermite_spline_t)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(curves::piecewise_t)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(curves::polynomial3_t)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(curves::exact_cubic3_t)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(curves::bezier3_t)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(curves::cubic_hermite_spline3_t)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(curves::piecewise3_t)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(curves::SO3Linear_t)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(curves::SE3Curve_t)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(curves::piecewise_SE3_t)


#endif  //_VARIABLES_PYTHON_BINDINGS
