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

variables_3_t fillWithZeros(const linear_variable_3_t& var, const std::size_t totalvar, const std::size_t i)
{
    variables_3_t res;
    std::vector<linear_variable_3_t>& vars = res.variables_;
    for (std::size_t idx = 0; idx < i; ++idx)
    {
        vars.push_back(linear_variable_3_t::Zero());
    }
    vars.push_back(var);
    for (std::size_t idx = i+1; idx < totalvar; ++idx)
    {
        vars.push_back(linear_variable_3_t::Zero());
    }
    return res;
}

std::vector<variables_3_t> computeLinearControlPoints(const point_list_t& matrices, const point_list_t& vectors)
{
    std::vector<variables_3_t> res;
    std::vector<linear_variable_3_t> variables = matrix3DFromEigenArray(matrices, vectors);
    // now need to fill all this with zeros...
    std::size_t totalvar = variables.size();
    for (std::size_t i = 0; i < totalvar; ++i)
    {
        res.push_back( fillWithZeros(variables[i],totalvar,i));
    }
    return res;
}

/*linear variable control points*/
bezier_linear_variable_t* wrapBezierLinearConstructor(const point_list_t& matrices, const point_list_t& vectors)
{
    std::vector<variables_3_t> asVector = computeLinearControlPoints(matrices, vectors);
    return new bezier_linear_variable_t(asVector.begin(), asVector.end(), 0., 1.) ;
}

bezier_linear_variable_t* wrapBezierLinearConstructorBounds(const point_list_t& matrices, const point_list_t& vectors, const real T_min, const real T_max)
{
    std::vector<variables_3_t> asVector = computeLinearControlPoints(matrices, vectors);
    return new bezier_linear_variable_t(asVector.begin(), asVector.end(), T_min, T_max) ;
}


LinearControlPointsHolder* wayPointsToLists(const bezier_linear_variable_t& self)
{
    typedef typename bezier_linear_variable_t::t_point_t t_point;
    typedef typename bezier_linear_variable_t::t_point_t::const_iterator cit_point;
    const t_point& wps = self.waypoints();
    // retrieve num variables.
    std::size_t dim = wps[0].variables_.size()*3;
    Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic> matrices (dim,wps.size() * 3);
    Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic> vectors  (dim,wps.size());
    int col = 0;
    for(cit_point cit = wps.begin(); cit != wps.end(); ++cit, ++col)
    {
        const std::vector<linear_variable_3_t>& variables = cit->variables_;
        int i = 0;
        for(std::vector<linear_variable_3_t>::const_iterator varit = variables.begin();
            varit != variables.end(); ++varit, i+=3)
        {
            vectors.block<3,1>(i,col)   =  varit->b_;
            matrices.block<3,3>(i,col*3) = varit->A_;
        }
    }
    LinearControlPointsHolder* res (new LinearControlPointsHolder);
    res->res = std::make_pair(matrices, vectors);
    return res;
}


// does not include end time
LinearBezierVector* split(const bezier_linear_variable_t& self,  const vectorX_t& times)
{
    LinearBezierVector* res (new LinearBezierVector);
    bezier_linear_variable_t current = self;
    real current_time = 0.;
    real tmp;
    for(int i = 0; i < times.rows(); ++i)
    {
        tmp =times[i];
        std::pair<bezier_linear_variable_t, bezier_linear_variable_t> pairsplit = current.split(tmp-current_time);
        res->beziers_.push_back(pairsplit.first);
        current = pairsplit.second;
        current_time += tmp-current_time;
    }
    res->beziers_.push_back(current);
    return res;
}

}
 // namespace curves
