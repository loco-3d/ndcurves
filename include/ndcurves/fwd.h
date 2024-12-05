/**
 * \file fwd.h
 * \brief forward declaration of all curves class
 * \author Pierre Fernbach
 * \version 0.1
 * \date 27/11/19
 *
 */

#ifndef CURVES_FWD_H
#define CURVES_FWD_H
#include <Eigen/Dense>
#include <memory>
#include <vector>

namespace ndcurves {

template <typename Time, typename Numeric, bool Safe, typename Point,
          typename Point_derivate>
struct curve_abc;

template <typename Time, typename Numeric, bool Safe, typename Point>
struct bezier_curve;

template <typename Time, typename Numeric, bool Safe, typename Point,
          typename Point_derivate>
struct constant_curve;

template <typename Time, typename Numeric, bool Safe, typename Point>
struct cubic_hermite_spline;

template <typename Time, typename Numeric, bool Safe, typename Point,
          typename T_Point, typename SplineBase>
struct exact_cubic;

template <typename Time, typename Numeric, bool Safe, typename Point,
          typename Point_derivate, typename CurveType>
struct piecewise_curve;

template <typename Time, typename Numeric, bool Safe, typename Point,
          typename T_Point>
struct polynomial;

template <typename Time, typename Numeric, bool Safe>
struct SE3Curve;

template <typename Time, typename Numeric, bool Safe, typename Point>
struct sinusoidal;

template <typename Time, typename Numeric, bool Safe>
struct SO3Linear;

template <typename Time, typename Numeric, bool Safe>
struct SO3Smooth;

template <typename Numeric>
struct Bern;

template <typename Point>
struct curve_constraints;

template <typename Numeric, bool Safe>
struct linear_variable;

template <typename Numeric>
struct quadratic_variable;

// typedef of the commonly used templates arguments :
// eigen types :
typedef Eigen::Matrix<double, 1, 1> point1_t;
typedef Eigen::Vector3d point3_t;
typedef Eigen::Matrix<double, 6, 1> point6_t;
typedef Eigen::VectorXd pointX_t;
typedef Eigen::Matrix<double, 3, 3> matrix3_t;
typedef Eigen::Matrix<double, 4, 4> matrix4_t;
typedef Eigen::Quaternion<double> quaternion_t;
typedef Eigen::Transform<double, 3, Eigen::Affine> transform_t;
typedef std::vector<point1_t, Eigen::aligned_allocator<point1_t>> t_point1_t;
typedef std::vector<point3_t, Eigen::aligned_allocator<point3_t>> t_point3_t;
typedef std::vector<pointX_t, Eigen::aligned_allocator<pointX_t>> t_pointX_t;
typedef Eigen::Ref<const matrix3_t> matrix3_t_cst_ref;

// abstract curves types:
typedef curve_abc<double, double, true, pointX_t, pointX_t>
    curve_abc_t;  // base abstract class
typedef curve_abc<double, double, true, point3_t, point3_t>
    curve_3_t;                          // generic class of curve of size 3
typedef curve_3_t curve_translation_t;  // generic class of a translation curve
typedef curve_abc<double, double, true, matrix3_t, point3_t>
    curve_rotation_t;  // templated class used for the rotation (return
                       // dimension are fixed)
typedef curve_abc<double, double, true, transform_t, point6_t>
    curve_SE3_t;  // templated abstract class used for all the se3 curves
                  // (return dimension are fixed)

// shared pointer to abstract types:
typedef std::shared_ptr<curve_abc_t> curve_ptr_t;
typedef std::shared_ptr<curve_3_t> curve3_ptr_t;
typedef std::shared_ptr<curve_rotation_t> curve_rotation_ptr_t;
typedef std::shared_ptr<curve_translation_t> curve_translation_ptr_t;
typedef std::shared_ptr<curve_SE3_t> curve_SE3_ptr_t;

// definition of all curves class with pointX as return type:
typedef polynomial<double, double, true, pointX_t, t_pointX_t> polynomial_t;
typedef exact_cubic<double, double, true, pointX_t, t_pointX_t, polynomial_t>
    exact_cubic_t;
typedef bezier_curve<double, double, true, pointX_t> bezier_t;
typedef linear_variable<double, true> linear_variable_t;
typedef bezier_curve<double, double, true, linear_variable_t>
    bezier_linear_variable_t;
typedef constant_curve<double, double, true, pointX_t, pointX_t> constant_t;
typedef cubic_hermite_spline<double, double, true, pointX_t>
    cubic_hermite_spline_t;
typedef piecewise_curve<double, double, true, pointX_t, pointX_t, curve_abc_t>
    piecewise_t;
typedef sinusoidal<double, double, true, pointX_t> sinusoidal_t;

// definition of all curves class with point3 as return type:
typedef polynomial<double, double, true, point3_t, t_point3_t> polynomial3_t;
typedef polynomial<double, double, true, point1_t, t_point1_t> polynomial1_t;
typedef exact_cubic<double, double, true, point3_t, t_point3_t, polynomial_t>
    exact_cubic3_t;
typedef bezier_curve<double, double, true, point3_t> bezier3_t;
typedef constant_curve<double, double, true, point3_t, point3_t> constant3_t;
typedef cubic_hermite_spline<double, double, true, point3_t>
    cubic_hermite_spline3_t;
typedef piecewise_curve<double, double, true, point3_t, point3_t, curve_3_t>
    piecewise3_t;

// special curves with return type fixed:
typedef SO3Smooth<double, double, true> SO3Smooth_t;
typedef SO3Linear<double, double, true> SO3Linear_t;
typedef SE3Curve<double, double, true> SE3Curve_t;
typedef piecewise_curve<double, double, true, transform_t, point6_t,
                        curve_SE3_t>
    piecewise_SE3_t;

}  // namespace ndcurves

#endif  // CURVES_FWD_H
