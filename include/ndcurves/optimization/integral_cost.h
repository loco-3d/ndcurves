/**
 * \file bezier_curve.h
 * \brief class allowing to create a Bezier curve of dimension 1 <= n <= 3.
 * \author Steve T.
 * \version 0.1
 * \date 06/17/2013
 */

#ifndef _CLASS_QUADRATIC_COST
#define _CLASS_QUADRATIC_COST

#include <Eigen/Core>

#include "ndcurves/optimization/definitions.h"
#include "ndcurves/optimization/details.h"

namespace ndcurves {
namespace optimization {

enum integral_cost_flag {
  DISTANCE = 0x000,
  VELOCITY = 0x001,
  ACCELERATION = 0x002,
  JERK = 0x003,
  FOURTH = 0x004,
  FIFTH = 0x005
};

template <typename Point, typename Numeric>
quadratic_variable<Numeric> compute_integral_cost_internal(
    const problem_data<Point, Numeric>& pData, const std::size_t num_derivate) {
  typedef bezier_curve<Numeric, Numeric, true, linear_variable<Numeric> >
      bezier_t;
  typedef typename bezier_t::t_point_t t_point_t;
  typedef typename t_point_t::const_iterator cit_point_t;
  bezier_t acc = pData.bezier->compute_derivate(num_derivate);
  const t_point_t& wps = acc.waypoints();
  quadratic_variable<Numeric> res(bezier_product<Point, Numeric, cit_point_t>(
      wps.begin(), wps.end(), wps.begin(), wps.end(), pData.dim_));
  return res;
}

template <typename Point, typename Numeric>
quadratic_variable<Numeric> compute_integral_cost(
    const problem_data<Point, Numeric>& pData, const integral_cost_flag flag) {
  std::size_t size = (std::size_t)(flag);
  return compute_integral_cost_internal<Point, Numeric>(pData, size);
}

}  // namespace optimization
}  // namespace ndcurves
#endif  //_CLASS_QUADRATIC_COST
