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

namespace curves {

template <typename Time, typename Numeric, bool Safe, typename Point , typename Point_derivate >
  struct curve_abc;

template <typename Time, typename Numeric, bool Safe,typename Point >
  struct bezier_curve;

template <typename Time, typename Numeric, bool Safe,typename Point >
  struct cubic_hermite_spline;

template <typename Time, typename Numeric, bool Safe, typename Point ,
            typename T_Point , typename SplineBase >
  struct exact_cubic;

template <typename Time, typename Numeric, bool Safe, typename Point,
            typename T_Point, typename Curve, typename Point_derivate>
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
}

#endif // CURVES_FWD_H
