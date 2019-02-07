#include "hpp/spline/bezier_curve.h"
#include "hpp/spline/linear_variable.h"

#include "python_definitions.h"

#include <vector>

#ifndef _VARIABLES_PYTHON_BINDINGS
#define _VARIABLES_PYTHON_BINDINGS


namespace spline
{
static const int dim = 3;
typedef spline::linear_variable<dim, real> linear_variable_3_t;
typedef spline::variables<linear_variable_3_t> variables_3_t;
typedef spline::bezier_curve  <real, real, dim, true, variables_3_t> bezier_linear_variable_t;


/*linear variable control points*/
bezier_linear_variable_t* wrapBezierLinearConstructor(const point_list_t& matrices, const point_list_t& vectors);

bezier_linear_variable_t* wrapBezierLinearConstructorBounds
    (const point_list_t& matrices, const point_list_t& vectors, const real ub);

typedef std::pair<Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic>,
                  Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic> > linear_points_t;

struct LinearControlPointsHolder
{
    linear_points_t res;
    Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic> A() {return res.first;}
    Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic> b() {return res.second;}
};

LinearControlPointsHolder* wayPointsToLists(const bezier_linear_variable_t& self);

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
LinearBezierVector* split(const bezier_linear_variable_t& self,  const vectorX_t& times);
} //namespace spline.


#endif //_VARIABLES_PYTHON_BINDINGS
