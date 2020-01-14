/**
 * \file bezier_curve.h
 * \brief class allowing to create a Bezier curve of dimension 1 <= n <= 3.
 * \author Steve T.
 * \version 0.1
 * \date 06/17/2013
 */

#ifndef _CLASS_LINEAR_PROBLEM_DETAILS
#define _CLASS_LINEAR_PROBLEM_DETAILS

#include <curves/bezier_curve.h>
#include <curves/linear_variable.h>
#include <curves/curve_constraint.h>
#include <curves/optimization/definitions.h>
#include <curves/bernstein.h>

namespace curves {
namespace optimization {
template <typename Point, typename Numeric, bool Safe = true>
struct problem_data {
  problem_data(const std::size_t dim) : bezier(0), dim_(dim) {}
  ~problem_data() {
    if (bezier) delete bezier;
  }

  typedef linear_variable<Numeric> var_t;
  typedef std::vector<var_t> T_var_t;
  typedef bezier_curve<Numeric, Numeric, true, linear_variable<Numeric> > bezier_t;

  std::vector<var_t> variables_;   // includes constant variables
  std::size_t numVariables;        // total number of variable (/ DIM for total size)
  std::size_t numControlPoints;    // total number of control Points (variables + waypoints) / DIM )
  std::size_t startVariableIndex;  // before that index, variables are constant
  std::size_t numStateConstraints;
  bezier_t* bezier;
  const std::size_t dim_;

  problem_data(const problem_data& other)
      : variables_(other.variables_),
        numVariables(other.numVariables),
        numControlPoints(other.numControlPoints),
        startVariableIndex(other.startVariableIndex),
        numStateConstraints(other.numStateConstraints),
        dim_(other.dim_) {
    const bezier_t& b = *other.bezier;
    bezier = new bezier_t(b.waypoints().begin(), b.waypoints().end(), b.T_min_, b.T_max_, b.mult_T_);
  }
};

inline std::size_t num_active_constraints(const constraint_flag& flag) {
  long lValue = (long)(flag);
  std::size_t iCount = 0;
  while (lValue != 0) {
    lValue = lValue & (lValue - 1);
    iCount++;
  }
  return (flag & NONE) ? iCount - 1 : iCount;
}

template <typename Numeric, typename LinearVar>
LinearVar fill_with_zeros(const LinearVar& var, const std::size_t i, const std::size_t startVariableIndex,
                          const std::size_t numVariables, const std::size_t Dim) {
  typedef Eigen::Matrix<Numeric, Eigen::Dynamic, Eigen::Dynamic> matrix_t;
  typename LinearVar::matrix_x_t B;
  B = matrix_t::Zero(Dim, numVariables * Dim);
  if (startVariableIndex <= i && i <= startVariableIndex + numVariables - 1 && var.size() > 0)
    B.block(0, Dim * (i - startVariableIndex), Dim, Dim) = var.B();
  return LinearVar(B, var.c());
}

template <typename Point, typename Numeric, typename Bezier, typename LinearVar>
Bezier* compute_linear_control_points(const problem_data<Point, Numeric>& pData,
                                      const std::vector<LinearVar>& linearVars, const Numeric totalTime) {
  std::vector<LinearVar> res;
  // now need to fill all this with zeros...
  std::size_t totalvar = linearVars.size();
  for (std::size_t i = 0; i < totalvar; ++i)
    res.push_back(fill_with_zeros<Numeric, LinearVar>(linearVars[i], i, pData.startVariableIndex, pData.numVariables,
                                                      pData.dim_));
  return new Bezier(res.begin(), res.end(), 0., totalTime);
}

template <typename Point, typename Numeric, bool Safe>
problem_data<Point, Numeric, Safe> setup_control_points(const problem_definition<Point, Numeric>& pDef) {
  typedef Numeric num_t;
  typedef Point point_t;
  typedef linear_variable<Numeric> var_t;
  typedef problem_data<Point, Numeric> problem_data_t;

  const std::size_t& degree = pDef.degree;
  const constraint_flag& flag = pDef.flag;

  const std::size_t numControlPoints = pDef.degree + 1;
  const std::size_t numActiveConstraints = num_active_constraints(flag);
  if (numActiveConstraints >= numControlPoints)
    throw std::runtime_error("In setup_control_points; too many constraints for the considered degree");

  problem_data_t problemData(pDef.dim_);
  typename problem_data_t::T_var_t& variables_ = problemData.variables_;

  std::size_t numConstants = 0;
  std::size_t i = 0;
  if (flag & INIT_POS) {
    variables_.push_back(var_t(pDef.init_pos));
    ++numConstants;
    ++i;
    if (flag & INIT_VEL) {
      point_t vel = pDef.init_pos + (pDef.init_vel / (num_t)degree) / pDef.totalTime;
      variables_.push_back(var_t(vel));
      ++numConstants;
      ++i;
      if (flag & INIT_ACC) {
        point_t acc = (pDef.init_acc / (num_t)(degree * (degree - 1))) / (pDef.totalTime * pDef.totalTime) + 2 * vel -
                      pDef.init_pos;
        ;
        variables_.push_back(var_t(acc));
        ++numConstants;
        ++i;
        if (flag & INIT_JERK) {
          point_t jerk = pDef.init_jerk * pDef.totalTime * pDef.totalTime * pDef.totalTime /
                             (num_t)(degree * (degree - 1) * (degree - 2)) +
                         3 * acc - 3 * vel + pDef.init_pos;
          variables_.push_back(var_t(jerk));
          ++numConstants;
          ++i;
        }
      }
    }
  }
  const std::size_t first_variable_idx = i;
  // variables
  for (; i + 4 < numControlPoints; ++i) variables_.push_back(var_t::Zero(pDef.dim_));
  // end constraints
  if (flag & END_POS) {
    if (flag & END_VEL) {
      point_t vel = pDef.end_pos - (pDef.end_vel / (num_t)degree) / pDef.totalTime;
      if (flag & END_ACC) {
        point_t acc = (pDef.end_acc / (num_t)(degree * (degree - 1))) / (pDef.totalTime) * (pDef.totalTime) + 2 * vel -
                      pDef.end_pos;
        if (flag & END_JERK) {
          point_t jerk = -pDef.end_jerk * pDef.totalTime * pDef.totalTime * pDef.totalTime /
                             (num_t)(degree * (degree - 1) * (degree - 2)) +
                         3 * acc - 3 * vel + pDef.end_pos;
          variables_.push_back(var_t(jerk));
          ++numConstants;
          ++i;
        } else
          while (i < numControlPoints - 3) {
            variables_.push_back(var_t::Zero(pDef.dim_));
            ++i;
          }
        variables_.push_back(var_t(acc));
        ++numConstants;
        ++i;
      } else
        while (i < numControlPoints - 2) {
          variables_.push_back(var_t::Zero(pDef.dim_));
          ++i;
        }
      variables_.push_back(var_t(vel));
      ++numConstants;
      ++i;
    } else {
      while (i < numControlPoints - 1) {
        variables_.push_back(var_t::Zero(pDef.dim_));
        ++i;
      }
    }
    variables_.push_back(var_t(pDef.end_pos));
    ++numConstants;
    ++i;
  }
  // add remaining variables (only if no end_pos constraints)
  for (; i < numControlPoints; ++i) variables_.push_back(var_t::Zero(pDef.dim_));

  if (numControlPoints < numConstants) {
    throw std::runtime_error("numControlPoints < numConstants");
  }
  if (numControlPoints != variables_.size()) {
    throw std::runtime_error("numControlPoints != variables_.size()");
  }

  problemData.numControlPoints = numControlPoints;
  problemData.numVariables = numControlPoints - numConstants;
  problemData.startVariableIndex = first_variable_idx;
  problemData.numStateConstraints = numActiveConstraints - problemData.numVariables;
  problemData.bezier =
      compute_linear_control_points<Point, Numeric, bezier_curve<Numeric, Numeric, true, var_t>, var_t>(
          problemData, variables_, pDef.totalTime);
  return problemData;
}

// TODO assumes constant are inside constraints...
template <typename Point, typename Numeric>
long compute_num_ineq_control_points(const problem_definition<Point, Numeric>& pDef,
                                     const problem_data<Point, Numeric>& pData) {
  typedef problem_definition<Point, Numeric> problem_definition_t;
  long rows(0);
  // rows depends on each constraint size, and the number of waypoints
  for (typename problem_definition_t::CIT_vector_x_t cit = pDef.inequalityVectors_.begin();
       cit != pDef.inequalityVectors_.end(); ++cit)
    rows += cit->rows() * pData.numControlPoints;
  return rows;
}

template <typename Point, typename Numeric>
std::vector<bezier_curve<Numeric, Numeric, true, linear_variable<Numeric> > > split(
    const problem_definition<Point, Numeric>& pDef, problem_data<Point, Numeric>& pData) {
  typedef linear_variable<Numeric> linear_variable_t;
  typedef bezier_curve<Numeric, Numeric, true, linear_variable_t> bezier_t;
  typedef std::vector<bezier_t> T_bezier_t;

  const Eigen::VectorXd& times = pDef.splitTimes_;
  T_bezier_t res;
  bezier_t& current = *pData.bezier;
  Numeric current_time = 0.;
  Numeric tmp;
  for (int i = 0; i < times.rows(); ++i) {
    tmp = times[i];
    std::pair<bezier_t, bezier_t> pairsplit = current.split(tmp - current_time);
    res.push_back(pairsplit.first);
    current = pairsplit.second;
    current_time += tmp - current_time;
  }
  res.push_back(current);
  return res;
}

template <typename Point, typename Numeric>
void initInequalityMatrix(const problem_definition<Point, Numeric>& pDef, problem_data<Point, Numeric>& pData,
                          quadratic_problem<Point, Numeric>& prob) {
  const std::size_t& Dim = pData.dim_;
  typedef problem_definition<Point, Numeric> problem_definition_t;
  typedef typename problem_definition_t::matrix_x_t matrix_x_t;
  typedef typename problem_definition_t::vector_x_t vector_x_t;
  typedef bezier_curve<Numeric, Numeric, true, linear_variable<Numeric> > bezier_t;
  typedef std::vector<bezier_t> T_bezier_t;
  typedef typename T_bezier_t::const_iterator CIT_bezier_t;
  typedef typename bezier_t::t_point_t t_point;
  typedef typename bezier_t::t_point_t::const_iterator cit_point;

  long cols = pData.numVariables * Dim;
  long rows = compute_num_ineq_control_points<Point, Numeric>(pDef, pData);
  prob.ineqMatrix = matrix_x_t::Zero(rows, cols);
  prob.ineqVector = vector_x_t::Zero(rows);

  if (pDef.inequalityMatrices_.size() == 0) return;

  // compute sub-bezier curves
  T_bezier_t beziers = split<Point, Numeric>(pDef, pData);

  if (pDef.inequalityMatrices_.size() != pDef.inequalityVectors_.size()) {
    throw std::invalid_argument("The sizes of the inequality matrices and vectors do not match.");
  }
  if (pDef.inequalityMatrices_.size() != beziers.size()) {
    throw std::invalid_argument("The sizes of the inequality matrices and the bezier degree do not match.");
  }

  long currentRowIdx = 0;
  typename problem_definition_t::CIT_matrix_x_t cmit = pDef.inequalityMatrices_.begin();
  typename problem_definition_t::CIT_vector_x_t cvit = pDef.inequalityVectors_.begin();
  // for each bezier split ..
  for (CIT_bezier_t bit = beziers.begin(); bit != beziers.end(); ++bit, ++cvit, ++cmit) {
    // compute vector of linear expressions of each control point
    const t_point& wps = bit->waypoints();
    // each control has a linear expression depending on all variables
    for (cit_point cit = wps.begin(); cit != wps.end(); ++cit) {
      prob.ineqMatrix.block(currentRowIdx, 0, cmit->rows(), cols) =
          (*cmit) * (cit->B());  // constraint inequality for current bezier * expression of control point
      prob.ineqVector.segment(currentRowIdx, cmit->rows()) = *cvit - (*cmit) * (cit->c());
      currentRowIdx += cmit->rows();
    }
  }
  assert(rows == currentRowIdx);  // we filled all the constraints - NB: leave assert for Debug tests
}

template <typename Point, typename Numeric, typename In>
quadratic_variable<Numeric> bezier_product(In PointsBegin1, In PointsEnd1, In PointsBegin2, In PointsEnd2,
                                           const std::size_t /*Dim*/) {
  typedef Eigen::Matrix<Numeric, Eigen::Dynamic, 1> vector_x_t;
  unsigned int nPoints1 = (unsigned int)(std::distance(PointsBegin1, PointsEnd1)),
               nPoints2 = (unsigned int)(std::distance(PointsBegin2, PointsEnd2));
  if (nPoints1 < 0 || nPoints2 < 0) {
    throw std::runtime_error("This should never happen because an unsigned int cannot go negative without underflowing.");
  }
  unsigned int deg1 = nPoints1 - 1, deg2 = nPoints2 - 1;
  unsigned int newDeg = (deg1 + deg2);
  // the integral of the primitive will simply be the last control points of the primitive,
  // divided by the degree of the primitive, newDeg. We will store this in matrices for bilinear terms,
  // and a vector for the linear terms, as well as another one for the constants.
  quadratic_variable<Numeric> res(vector_x_t::Zero(PointsBegin1->B().cols()));
  // depending on the index, the fraction coefficient of the bernstein polynom
  // is either the fraction given by  (i+j)/ (deg1+deg2), or 1 - (i+j)/ (deg1+deg2).
  // The trick is that the condition is given by whether the current index in
  // the combinatorial is odd or even.
  // time parametrization is not relevant for the cost

  Numeric ratio;
  for (unsigned int i = 0; i < newDeg + 1; ++i) {
    unsigned int j = i > deg2 ? i - deg2 : 0;
    for (; j < std::min(deg1, i) + 1; ++j) {
      ratio = (Numeric)(bin(deg1, j) * bin(deg2, i - j)) / (Numeric)(bin(newDeg, i));
      In itj = PointsBegin1 + j;
      In iti = PointsBegin2 + (i - j);
      res += ((*itj) * (*iti)) * ratio;
    }
  }
  return res / (newDeg + 1);
}

inline constraint_flag operator~(constraint_flag a) {
  return static_cast<constraint_flag>(~static_cast<const int>(a));
}

inline constraint_flag operator|(constraint_flag a, constraint_flag b) {
  return static_cast<constraint_flag>(static_cast<const int>(a) | static_cast<const int>(b));
}

inline constraint_flag operator&(constraint_flag a, constraint_flag b) {
  return static_cast<constraint_flag>(static_cast<const int>(a) & static_cast<const int>(b));
}

inline constraint_flag operator^(constraint_flag a, constraint_flag b) {
  return static_cast<constraint_flag>(static_cast<const int>(a) ^ static_cast<const int>(b));
}

inline constraint_flag& operator|=(constraint_flag& a, constraint_flag b) {
  return (constraint_flag&)((int&)(a) |= static_cast<const int>(b));
}

inline constraint_flag& operator&=(constraint_flag& a, constraint_flag b) {
  return (constraint_flag&)((int&)(a) &= static_cast<const int>(b));
}

inline constraint_flag& operator^=(constraint_flag& a, constraint_flag b) {
  return (constraint_flag&)((int&)(a) ^= static_cast<const int>(b));
}

}  // namespace optimization
}  // namespace curves
#endif  //_CLASS_LINEAR_PROBLEM_DETAILS
