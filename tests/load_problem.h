
#ifndef _CLASS_LOAD_TEST_PROBLEMS
#define _CLASS_LOAD_TEST_PROBLEMS


#include "curves/exact_cubic.h"
#include "curves/bezier_curve.h"
#include "curves/helpers/effector_spline.h"
#include "curves/helpers/effector_spline_rotation.h"
#include "curves/optimization/quadratic_problem.h"
#include "curves/optimization/integral_cost.h"
#include "curves/optimization/details.h"

#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>

namespace curves {

typedef Eigen::Vector3d point_t;
typedef std::vector<point_t,Eigen::aligned_allocator<point_t> >  t_point_t;
typedef Eigen::VectorXd pointX_t;
typedef std::pair<double, pointX_t> Waypoint;
typedef std::vector<Waypoint> T_Waypoint;


typedef Eigen::Matrix<double,1,1> point_one;
typedef std::pair<double, point_one> WaypointOne;
typedef std::vector<WaypointOne> T_WaypointOne;

namespace optimization
{
typedef curve_constraints<point_t> constraint_linear;
typedef linear_variable<double,true> linear_variable_t;
typedef std::vector<linear_variable_t> T_linear_variable_t;
typedef T_linear_variable_t::const_iterator CIT_linear_variable_t;
typedef std::pair<std::size_t, std::size_t >   pair_size_t;
typedef std::pair<T_linear_variable_t, pair_size_t > var_pair_t;
typedef problem_data<point_t, double> problem_data_t;
typedef problem_definition<point_t, double> problem_definition_t;
typedef quadratic_problem<point_t, double> problem_t;


#define MAXBUFSIZE  ((int) 1e6)

Eigen::MatrixXd readMatrix(std::ifstream& infile)
{
    int cols = 0, rows = 0;
    double buff[MAXBUFSIZE];

    // Read numbers from file into buffer.
    //ifstream infile;
    //infile.open(filename);
    std::string line = "noise";
    while (!infile.eof() && !line.empty())
    {
        std::getline(infile, line);

        int temp_cols = 0;
        std::stringstream stream(line);
        while(! stream.eof())
            stream >> buff[cols*rows+temp_cols++];

        if (temp_cols == 0)
            continue;

        if (cols == 0)
            cols = temp_cols;

        rows++;
    }
    //infile.close();
    rows--;

    // Populate matrix with numbers.
    Eigen::MatrixXd result(rows,cols);
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            result(i,j) = buff[ cols*i+j ];

    return result;
}

problem_definition_t loadproblem(const std::string& filename)
{
    problem_definition_t pDef(3);
    std::ifstream in (filename.c_str());
    if (!in.is_open())
        throw std::runtime_error("cant open filename");
    //first line is degree totaltime flag
    Eigen::Vector3d degTimeFlag = readMatrix(in);
    pDef.degree = (std::size_t)(degTimeFlag[0]);
    pDef.totalTime =degTimeFlag[1];
    pDef.flag = (constraint_flag)(degTimeFlag[2]);
    //Then startpos then empty line
    pDef.init_pos = readMatrix(in);
    //Then endpos then empty line
    pDef.end_pos = readMatrix(in);
    //Then splittimes then empty line
    pDef.splitTimes_ = readMatrix(in);
    // The inequality matrices, empty line, inequality vector as many times
    for (int i = 0; i< pDef.splitTimes_.rows()+1; ++i)
    {
        pDef.inequalityMatrices_.push_back(readMatrix(in));
        pDef.inequalityVectors_.push_back(readMatrix(in));
    }
    in.close();
    return pDef;
    // TODO curve constraints
}

}
}


#endif //_CLASS_LOAD_TEST_PROBLEMS
