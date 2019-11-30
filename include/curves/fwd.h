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
#include <vector>
#include <boost/smart_ptr/shared_ptr.hpp>

namespace curves {

template <typename Time, typename Numeric, bool Safe, typename Point , typename Point_derivate >
  struct curve_abc;

template <typename Time, typename Numeric, bool Safe,typename Point >
  struct bezier_curve;

template <typename Time, typename Numeric, bool Safe,typename Point >
  struct cubic_hermite_spline;

template <typename Time, typename Numeric, bool Safe, typename Point,
          typename T_Point,typename SplineBase >
  struct exact_cubic;

template <typename Time, typename Numeric, bool Safe, typename Point,
            typename Point_derivate>
  struct piecewise_curve;

template <typename Time, typename Numeric, bool Safe,typename Point, typename T_Point>
  struct polynomial;

template <typename Time, typename Numeric, bool Safe>
  struct SE3Curve;

template <typename Time, typename Numeric, bool Safe>
  struct SO3Linear;

template <typename Numeric>
  struct Bern;

template <typename Point>
  struct curve_constraints;

template <typename Numeric, bool Safe>
  struct linear_variable;

template <typename Numeric>
  struct quadratic_variable;

// typedef of the commonly used templates arguments :

typedef Eigen::Vector3d point3_t;
typedef Eigen::Matrix<double, 6, 1> point6_t;
typedef Eigen::VectorXd pointX_t;
typedef Eigen::Matrix<double, 3, 3> matrix3_t;
typedef Eigen::Quaternion<double> quaternion_t;
typedef Eigen::Transform<double, 3, Eigen::Affine> transform_t;
typedef std::vector<pointX_t, Eigen::aligned_allocator<pointX_t> > t_pointX_t;
typedef curve_abc<double, double, true, pointX_t,pointX_t> curve_abc_t; // base abstract class
typedef curve_abc<double, double, true, matrix3_t, point3_t> curve_rotation_t;  // templated class used for the rotation (return dimension are fixed)
typedef boost::shared_ptr<curve_abc_t> curve_ptr_t;
typedef boost::shared_ptr<curve_rotation_t> curve_rotation_ptr_t;
typedef polynomial<double, double, true, pointX_t, t_pointX_t> polynomial_t;
typedef exact_cubic<double, double, true, pointX_t,
  std::vector<pointX_t, Eigen::aligned_allocator<pointX_t> >, polynomial_t> exact_cubic_t;
typedef bezier_curve<double, double, true, pointX_t> bezier_curve_t;
typedef SO3Linear<double, double, true> SO3Linear_t;
typedef SE3Curve<double, double, true> SE3Curve_t;
typedef curve_abc<double, double, true, transform_t, point6_t> curve_SE3_t;  // templated abstract class used for all the se3 curves (return dimension are fixed)
typedef boost::shared_ptr<curve_SE3_t> curve_SE3_ptr_t;
typedef cubic_hermite_spline<double, double, true, pointX_t> cubic_hermite_spline_t;
typedef piecewise_curve <double, double, true, pointX_t,pointX_t> piecewise_curve_t;
typedef piecewise_curve <double, double, true,  transform_t, point6_t> piecewise_SE3_curve_t;

}

#endif // CURVES_FWD_H
