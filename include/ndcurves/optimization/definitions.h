/**
 * \file definitions.h
 * \brief utils for defining optimization problems
 * \author Steve T.
 * \version 0.1
 * \date 06/17/2013
 */

#ifndef _CLASS_DEFINITIONS_H
#define _CLASS_DEFINITIONS_H

#include <ndcurves/bezier_curve.h>
#include <ndcurves/curve_constraint.h>
#include <ndcurves/linear_variable.h>
#include <ndcurves/quadratic_variable.h>

namespace ndcurves {
namespace optimization {

enum constraint_flag {
  INIT_POS = 0x001,
  INIT_VEL = 0x002,
  INIT_ACC = 0x004,
  INIT_JERK = 0x008,
  END_POS = 0x010,
  END_VEL = 0x020,
  END_ACC = 0x040,
  END_JERK = 0x080,
  ALL = 0x0ff,
  NONE = 0x100
};

template <typename Point, typename Numeric>
struct quadratic_problem {
  Eigen::Matrix<Numeric, Eigen::Dynamic, Eigen::Dynamic> ineqMatrix;
  Eigen::Matrix<Numeric, Eigen::Dynamic, 1> ineqVector;
  quadratic_variable<Numeric> cost;
};

template <typename Point, typename Numeric>
struct problem_definition : public curve_constraints<Point> {
  typedef Point point_t;
  typedef Numeric num_t;
  typedef curve_constraints<point_t> curve_constraints_t;
  typedef Eigen::Matrix<num_t, Eigen::Dynamic, 1> vector_x_t;
  typedef Eigen::Matrix<num_t, Eigen::Dynamic, Eigen::Dynamic> matrix_x_t;
  typedef std::vector<matrix_x_t, Eigen::aligned_allocator<matrix_x_t> >
      T_matrix_x_t;
  typedef std::vector<vector_x_t, Eigen::aligned_allocator<vector_x_t> >
      T_vector_x_t;
  typedef typename T_matrix_x_t::const_iterator CIT_matrix_x_t;
  typedef typename T_vector_x_t::const_iterator CIT_vector_x_t;

  problem_definition(const std::size_t dim)
      : curve_constraints_t(dim),
        flag(NONE),
        init_pos(point_t::Zero(dim)),
        end_pos(point_t::Zero(dim)),
        degree(5),
        totalTime(1.),
        splitTimes_(vector_x_t::Zero(0)),
        dim_(dim) {}

  problem_definition(const curve_constraints_t& parent)
      : curve_constraints_t(parent),
        flag(NONE),
        init_pos(point_t::Zero(parent.dim_)),
        end_pos(point_t::Zero(parent.dim_)),
        degree(5),
        totalTime(1.),
        splitTimes_(vector_x_t::Zero(0)),
        dim_(parent.dim_) {}

  constraint_flag flag;
  point_t init_pos;
  point_t end_pos;
  std::size_t degree;
  num_t totalTime;
  vector_x_t splitTimes_;
  T_matrix_x_t inequalityMatrices_;  // must be of size (splitTimes_ + 1)
  T_vector_x_t inequalityVectors_;   // must be of size (splitTimes_ + 1)
  const std::size_t dim_;
};

}  // namespace optimization
}  // namespace ndcurves
#endif  //_CLASS_DEFINITIONS_H
