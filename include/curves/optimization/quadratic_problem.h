/**
* \file bezier_curve.h
* \brief class allowing to create a Bezier curve of dimension 1 <= n <= 3.
* \author Steve T.
* \version 0.1
* \date 06/17/2013
*/


#ifndef _CLASS_LINEAR_PROBLEM
#define _CLASS_LINEAR_PROBLEM

#include "curves/optimization/definitions.h"
#include "curves/optimization/details.h"
#include "curves/optimization/integral_cost.h"

#include <Eigen/Core>

namespace curves
{
namespace  optimization
{

template<typename Point, int Dim, typename Numeric>
quadratic_problem<Point, Dim, Numeric> generate_problem
    (const problem_definition<Point, Dim, Numeric>& pDef, const quadratic_variable<Numeric>& cost)
{
    quadratic_problem<Point, Dim, Numeric> prob;
    problem_data<Point, Dim, Numeric> pData = setup_control_points<Point, Dim, Numeric>(pDef);
    initInequalityMatrix<Point, Dim, Numeric>(pDef,pData,prob);
    prob.cost = cost;
    return prob;
}

template<typename Point, int Dim, typename Numeric>
quadratic_problem<Point, Dim, Numeric> generate_problem
    (const problem_definition<Point, Dim, Numeric>& pDef, const integral_cost_flag costFlag)
{
    quadratic_problem<Point, Dim, Numeric> prob;
    problem_data<Point, Dim, Numeric> pData = setup_control_points<Point, Dim, Numeric>(pDef);
    initInequalityMatrix<Point, Dim, Numeric>(pDef,pData,prob);
    prob.cost = compute_integral_cost<Point, Dim, Numeric>(pData, costFlag);
    return prob;
}
} // namespace optimization
} // namespace curves
#endif //_CLASS_LINEAR_PROBLEM

