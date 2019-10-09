#include "python_variables.h"
#include "python_definitions.h"

#include <Eigen/Core>

namespace curves
{
  std::vector<linear_variable_3_t> matrix3DFromEigenArray(const point_list_t& matrices, const point_list_t& vectors)
  {
    assert(vectors.cols() * 3  == matrices.cols() ) ;
    std::vector<linear_variable_3_t> res;
    for(int i =0;i<vectors.cols();++i)
    {
      res.push_back(linear_variable_3_t(matrices.block<3,3>(0,i*3), vectors.col(i)));
    }
    return res;
  }

  linear_variable_3_t fillWithZeros(const linear_variable_3_t& var, const std::size_t totalvar, const std::size_t i)
  {
      linear_variable_3_t::matrix_x_t B(linear_variable_3_t::matrix_x_t::Zero(dim,totalvar*dim));
      B.block(0,dim*i,dim,dim) = var.B();
      return linear_variable_3_t (B,var.c());
  }

  std::vector<linear_variable_3_t> computeLinearControlPoints(const point_list_t& matrices, const point_list_t& vectors)
  {
      std::vector<linear_variable_3_t> res;
      std::vector<linear_variable_3_t> variables = matrix3DFromEigenArray(matrices, vectors);
      // now need to fill all this with zeros...
      std::size_t totalvar = variables.size();
      for (std::size_t i = 0; i < totalvar; ++i)
          res.push_back( fillWithZeros(variables[i],totalvar,i));
      return res;
  }

  /*linear variable control points*/
  bezier_linear_variable_t* wrapBezierLinearConstructor(const point_list_t& matrices, const point_list_t& vectors)
  {
      std::vector<linear_variable_3_t> asVector = computeLinearControlPoints(matrices, vectors);
      return new bezier_linear_variable_t(asVector.begin(), asVector.end()) ;
  }

  bezier_linear_variable_t* wrapBezierLinearConstructorBounds(const point_list_t& matrices, const point_list_t& vectors, const real T_min, const real T_max)
  {
      std::vector<linear_variable_3_t> asVector = computeLinearControlPoints(matrices, vectors);
      return new bezier_linear_variable_t(asVector.begin(), asVector.end(), T_min, T_max) ;
  }


  Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic> cost_t_quad(const quadratic_variable_t& p)
  {
      Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic> A = p.A();
      return A;
  }
  Eigen::Matrix<real, Eigen::Dynamic, 1> cost_t_linear(const quadratic_variable_t & p)
  {
      Eigen::Matrix<real, Eigen::Dynamic, 1> b = p.b();
      return b;
  }
  real cost_t_constant(const quadratic_variable_t & p)
  {
      return p.c();
  }


  matrix_pair*
          wayPointsToLists(const bezier_linear_variable_t& self)
  {
      typedef typename bezier_linear_variable_t::t_point_t t_point;
      typedef typename bezier_linear_variable_t::t_point_t::const_iterator cit_point;
      const t_point& wps = self.waypoints();
      // retrieve num variables.
      std::size_t dim = wps[0].B().cols();
      Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic> matrices (dim,wps.size() * 3);
      Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic> vectors =
              Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic>::Zero(3,wps.size());
      int i = 0;
      for(cit_point cit = wps.begin(); cit != wps.end(); ++cit, ++i)
      {
          matrices.block(0,i*3,dim,3) = cit->B().transpose();
          vectors.block<3,1>(0,i)   =  cit->c();
      }
      matrix_pair* res (new matrix_pair(matrices, vectors));
      //res->res = std::make_pair(matrices, vectors);
      return res;
  }

  // does not include end time
  LinearBezierVector* split_py(const bezier_linear_variable_t& self,  const vectorX_t& times)
  {
      LinearBezierVector* res (new LinearBezierVector);
      bezier_linear_variable_t current = self;
      for(int i = 0; i < times.rows(); ++i)
      {
          std::pair<bezier_linear_variable_t, bezier_linear_variable_t> pairsplit = current.split(times[i]);
          res->beziers_.push_back(pairsplit.first);
          current = pairsplit.second;
      }
      res->beziers_.push_back(current);
      return res;
  }
} // namespace curves
