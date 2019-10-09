/**
* \file definitions.h
* \brief utils for defining optimization problems
* \author Steve T.
* \version 0.1
* \date 06/17/2013
*/


#ifndef _CLASS_DEFINITIONS_H
#define _CLASS_DEFINITIONS_H


#include <curves/bezier_curve.h>
#include <curves/linear_variable.h>
#include <curves/quadratic_variable.h>
#include <curves/curve_constraint.h>


namespace curves
{
namespace  optimization
{

enum constraint_flag{
    INIT_POS  = 0x001,
    INIT_VEL  = 0x002,
    INIT_ACC  = 0x004,
    INIT_JERK = 0x008,
    END_POS   = 0x010,
    END_VEL   = 0x020,
    END_ACC   = 0x040,
    END_JERK  = 0x080,
    ALL       = 0x0ff,
    NONE      = 0x100
};



template<typename Point, typename Numeric>
struct quadratic_problem
{
    Eigen::Matrix<Numeric,Eigen::Dynamic, Eigen::Dynamic> ineqMatrix;
    Eigen::Matrix<Numeric,Eigen::Dynamic, 1> ineqVector;
    quadratic_variable<Numeric> cost;
};


template<typename Point, typename Numeric>
struct problem_definition
{
    typedef Point  point_t;
    typedef Numeric  num_t;
    typedef curve_constraints<point_t> curve_constraints_t;
    typedef Eigen::Matrix< num_t , Eigen::Dynamic , 1> vector_x_t;
    typedef Eigen::Matrix< num_t , Eigen::Dynamic , Eigen::Dynamic> matrix_x_t;
    typedef Eigen::Matrix< num_t , Eigen::Dynamic , Eigen::Dynamic> matrix_dim_t;
    typedef std::vector<matrix_dim_t, Eigen::aligned_allocator<matrix_dim_t> > T_matrix_dim_t;
    typedef std::vector<vector_x_t, Eigen::aligned_allocator<vector_x_t> > T_vector_x_t;
    typedef typename T_matrix_dim_t::const_iterator CIT_matrix_dim_t;
    typedef typename T_vector_x_t::const_iterator CIT_vector_x_t;

    problem_definition(const int dim)
        : flag(NONE)
        , start(point_t::Zero(dim))
        , end(point_t::Zero(dim))
        , curveConstraints(dim)
        , degree(5)
        , totalTime(1.)
        , splitTimes_(vector_x_t::Zero(0))
        , dim_(dim){}


    constraint_flag flag;
    point_t start;
    point_t end;
    curve_constraints_t curveConstraints;
    std::size_t degree;
    num_t totalTime;
    vector_x_t splitTimes_;
    T_matrix_dim_t inequalityMatrices_; // must be of size (splitTimes_ + 1)
    T_vector_x_t    inequalityVectors_;  // must be of size (splitTimes_ + 1)
    const int dim_;
};

} // namespace optimization
} // namespace curves
#endif //_CLASS_DEFINITIONS_H

