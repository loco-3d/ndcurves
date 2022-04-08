/**
 * \file bezier_curve.h
 * \brief class allowing to create a Bezier curve of dimension 1 <= n <= 3.
 * \author Steve T.
 * \version 0.1
 * \date 06/17/2013
 */

#ifndef _CLASS_LINEAR_PROBLEM
#define _CLASS_LINEAR_PROBLEM

#include <Eigen/Core>

#include "ndcurves/optimization/definitions.h"
#include "ndcurves/optimization/details.h"
#include "ndcurves/optimization/integral_cost.h"

namespace ndcurves {
namespace optimization {

template <typename Point, typename Numeric, bool Safe>
quadratic_problem<Point, Numeric> generate_problem(
    const problem_definition<Point, Numeric>& pDef,
    const quadratic_variable<Numeric>& cost) {
  quadratic_problem<Point, Numeric> prob;
  problem_data<Point, Numeric> pData =
      setup_control_points<Point, Numeric, Safe>(pDef);
  initInequalityMatrix<Point, Numeric>(pDef, pData, prob);
  prob.cost = cost;
  return prob;
}

template <typename Point, typename Numeric, bool Safe>
quadratic_problem<Point, Numeric> generate_problem(
    const problem_definition<Point, Numeric>& pDef,
    const integral_cost_flag costFlag) {
  quadratic_problem<Point, Numeric> prob;
  problem_data<Point, Numeric> pData =
      setup_control_points<Point, Numeric, Safe>(pDef);
  initInequalityMatrix<Point, Numeric>(pDef, pData, prob);
  prob.cost = compute_integral_cost<Point, Numeric>(pData, costFlag);
  return prob;
}
}  // namespace optimization
}  // namespace ndcurves
#endif  //_CLASS_LINEAR_PROBLEM
