#include "curves/bezier_curve.h"
#include "curves/linear_variable.h"

#include "python_definitions.h"

#include <vector>

#ifndef _VARIABLES_PYTHON_BINDINGS
#define _VARIABLES_PYTHON_BINDINGS


namespace curves
{
  static const int dim = 3;
  typedef linear_variable<dim, real> linear_variable_3_t;
  typedef quadratic_variable<real> quadratic_variable_t;
  typedef bezier_curve  <real, real, true, linear_variable_3_t> bezier_linear_variable_t;

  /*linear variable control points*/
  bezier_linear_variable_t* wrapBezierLinearConstructor(const point_list_t& matrices, const point_list_t& vectors);
  bezier_linear_variable_t* wrapBezierLinearConstructorBounds
  (const point_list_t& matrices, const point_list_t& vectors, const real T_min, const real T_max);

  typedef std::pair<Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic>,
  Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic> > linear_points_t;

  struct matrix_pair
  {
      linear_points_t res;
      Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic> A() {return res.first;}
      Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic> b() {return res.second;}
  };


  Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic> cost_t_quad(const quadratic_variable_t& p);
  Eigen::Matrix<real, Eigen::Dynamic, 1> cost_t_linear(const quadratic_variable_t & p);
  real cost_t_constant(const quadratic_variable_t & p);

  matrix_pair* wayPointsToLists(const bezier_linear_variable_t& self);

  struct LinearBezierVector
  {
    std::vector<bezier_linear_variable_t> beziers_;
    std::size_t size() {return beziers_.size();}
    bezier_linear_variable_t* at(std::size_t i)
    {
      assert (i<size());
      return new bezier_linear_variable_t(beziers_[i]);
    }
  };

  // does not include end time
  LinearBezierVector* split_py(const bezier_linear_variable_t& self,  const vectorX_t& times);
} //namespace curve.


#endif //_VARIABLES_PYTHON_BINDINGS
