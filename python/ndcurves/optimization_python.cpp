#include <boost/python.hpp>
#include <boost/python/bases.hpp>
#include <boost/python/enum.hpp>

#include "archive_python_binding.h"
#include "namespace.h"
#include "ndcurves/optimization/definitions.h"
#include "ndcurves/optimization/quadratic_problem.h"
#include "python_variables.h"

namespace ndcurves {
namespace optimization {
namespace python {
static const bool safe = true;
typedef problem_definition<pointX_t, real> problem_definition_t;
typedef problem_data<pointX_t, real> problem_data_t;
typedef quadratic_problem<pointX_t, real> quadratic_problem_t;

problem_data_t setup_control_points_t(problem_definition_t& pDef) {
  problem_data_t pData = setup_control_points<pointX_t, real, safe>(pDef);
  return pData;  // return new problem_data_t(pData);
}

quadratic_variable_t problem_t_cost(const quadratic_problem_t& p) {
  return p.cost;
}
Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic> problem_t_ineqMatrix(
    const quadratic_problem_t& p) {
  return p.ineqMatrix;
}
Eigen::Matrix<real, Eigen::Dynamic, 1> problem_t_ineqVector(
    const quadratic_problem_t& p) {
  return p.ineqVector;
}

Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic> cost_t_quad(
    const quadratic_variable_t& p) {
  Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic> A = p.A();
  return A;
}
Eigen::Matrix<real, Eigen::Dynamic, 1> cost_t_linear(
    const quadratic_variable_t& p) {
  Eigen::Matrix<real, Eigen::Dynamic, 1> b = p.b();
  return b;
}
real cost_t_constant(const quadratic_variable_t& p) { return p.c(); }

quadratic_problem_t generate_problem_t(const problem_definition_t& pDef,
                                       const quadratic_variable_t& c) {
  return generate_problem<pointX_t, real, true>(pDef, c);
}

quadratic_problem_t generate_integral_problem_t(
    const problem_definition_t& pDef, const integral_cost_flag c) {
  return generate_problem<problem_definition_t::point_t, real, true>(pDef, c);
}

void set_pd_flag(problem_definition_t* pDef, const int flag) {
  pDef->flag = (constraint_flag)(flag);
}
void set_start(problem_definition_t* pDef, const pointX_t& val) {
  pDef->init_pos = val;
}
void set_end(problem_definition_t* pDef, const pointX_t& val) {
  pDef->end_pos = val;
}
void set_degree(problem_definition_t* pDef, const std::size_t val) {
  pDef->degree = val;
}
void set_total_time(problem_definition_t* pDef, const double val) {
  pDef->totalTime = val;
}
void set_split_time(problem_definition_t* pDef, const Eigen::VectorXd& val) {
  pDef->splitTimes_ = val;
}
Eigen::VectorXd get_split_times(const problem_definition_t* pDef) {
  return pDef->splitTimes_;
}

constraint_flag get_pd_flag(const problem_definition_t* pDef) {
  return pDef->flag;
}
Eigen::VectorXd get_start(const problem_definition_t* pDef) {
  return pDef->init_pos;
}
Eigen::VectorXd get_end(const problem_definition_t* pDef) {
  return pDef->end_pos;
}
std::size_t get_degree(const problem_definition_t* pDef) {
  return pDef->degree;
}
double get_total_time(const problem_definition_t* pDef) {
  return pDef->totalTime;
}

matrix_pair* get_ineq_at(const problem_definition_t* pDef,
                         const std::size_t idx) {
  if (idx > pDef->inequalityMatrices_.size() - 1)
    throw std::runtime_error(
        "required id is beyond number of inequality matrices");
  matrix_pair* res = new matrix_pair(pDef->inequalityMatrices_[idx],
                                     pDef->inequalityVectors_[idx]);
  return res;
}
bool del_ineq_at(problem_definition_t* pDef, const std::size_t idx) {
  if (idx > pDef->inequalityMatrices_.size() - 1) return false;
  pDef->inequalityMatrices_.erase(pDef->inequalityMatrices_.begin() + idx - 1);
  pDef->inequalityVectors_.erase(pDef->inequalityVectors_.begin() + idx - 1);
  return true;
}
bool add_ineq_at(problem_definition_t* pDef, const Eigen::MatrixXd ineq,
                 const Eigen::VectorXd vec) {
  if (ineq.rows() != vec.rows())
    throw std::runtime_error(
        "ineq vector and matrix do not have the same number of rows");
  if (!(pDef->inequalityMatrices_.empty()) &&
      ineq.cols() != pDef->inequalityMatrices_.back().cols())
    throw std::runtime_error(
        "inequality matrix does not have the same variable dimension as "
        "existing matrices");
  pDef->inequalityMatrices_.push_back(ineq);
  pDef->inequalityVectors_.push_back(vec);
  return true;
}

bezier_linear_variable_t* pDataBezier(const problem_data_t* pData) {
  const bezier_linear_variable_t& b = *pData->bezier;
  return new bezier_linear_variable_t(
      b.waypoints().begin(), b.waypoints().end(), b.min(), b.max(), b.mult_T_);
}

problem_definition_t* wrapProblemDefinitionConstructor(
    const curve_constraints_t* c) {
  return new problem_definition_t(*c);
}

void exposeOptimization() {
  // using the optimization scope
  bp::scope current_scope =
      ndcurves::python::getOrCreatePythonNamespace("optimization");
  /** BEGIN enums**/
  bp::enum_<constraint_flag>("constraint_flag")
      .value("INIT_POS", INIT_POS)
      .value("INIT_VEL", INIT_VEL)
      .value("INIT_ACC", INIT_ACC)
      .value("INIT_JERK", INIT_JERK)
      .value("END_POS", END_POS)
      .value("END_VEL", END_VEL)
      .value("END_ACC", END_ACC)
      .value("END_JERK", END_JERK)
      .value("ALL", ALL)
      .value("NONE", NONE)
      .export_values();

  bp::enum_<integral_cost_flag>("integral_cost_flag")
      .value("DISTANCE", DISTANCE)
      .value("VELOCITY", VELOCITY)
      .value("ACCELERATION", ACCELERATION)
      .value("JERK", JERK)
      .value("FOURTH", FOURTH)
      .value("FIFTH", FIFTH)
      .export_values();
  /** END enum**/

  bp::class_<quadratic_problem_t>("quadratic_problem", bp::init<>())
      .add_property("cost", &problem_t_cost)
      .add_property("A", &problem_t_ineqMatrix)
      .add_property("b", &problem_t_ineqVector);

  bp::def("setup_control_points", &setup_control_points_t);
  bp::def("generate_problem", &generate_problem_t);
  bp::def("generate_integral_problem", &generate_integral_problem_t);

  bp::class_<problem_data_t>("problem_data", bp::no_init)
      .def("bezier", &pDataBezier,
           bp::return_value_policy<bp::manage_new_object>())
      .def_readonly("numControlPoints", &problem_data_t::numControlPoints)
      .def_readonly("numVariables", &problem_data_t::numVariables)
      .def_readonly("startVariableIndex", &problem_data_t::startVariableIndex)
      .def_readonly("numStateConstraints",
                    &problem_data_t::numStateConstraints);

  bp::class_<problem_definition_t, bp::bases<curve_constraints_t> >(
      "problem_definition", bp::init<int>())
      .def("__init__", bp::make_constructor(&wrapProblemDefinitionConstructor))
      .add_property("flag", &get_pd_flag, &set_pd_flag)
      .add_property("init_pos", &get_start, &set_start)
      .add_property("end_pos", &get_end, &set_end)
      .add_property("degree", &get_degree, &set_degree)
      .add_property("totalTime", &get_total_time, &set_total_time)
      .add_property("splits", &get_split_times, &set_split_time)
      .def("inequality", &get_ineq_at,
           bp::return_value_policy<bp::manage_new_object>())
      .def("removeInequality", &del_ineq_at)
      .def("addInequality", &add_ineq_at);
}

}  // namespace python
}  // namespace optimization
}  // namespace ndcurves
