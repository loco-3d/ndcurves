#include "spline/bezier_curve.h"

#include <vector>

#include <eigenpy/memory.hpp>
#include <eigenpy/eigenpy.hpp>

#include <boost/python.hpp>

/*** TEMPLATE SPECIALIZATION FOR PYTHON ****/
typedef double real;
typedef Eigen::Vector3d point_t;
typedef Eigen::Matrix<real, Eigen::Dynamic,3> point_list_t;
typedef std::vector<point_t,Eigen::aligned_allocator<point_t> >  t_point_t;
typedef std::pair<real, point_t> Waypoint;
typedef std::vector<Waypoint> T_Waypoint;

typedef spline::bezier_curve  <real, real, 3, true, point_t> bezier_t;
/*** TEMPLATE SPECIALIZATION FOR PYTHON ****/

EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(bezier_t)

namespace spline
{
using namespace boost::python;

t_point_t vectorFromEigenArray(const point_list_t& array)
{
    t_point_t res;
    for(int i =0;i<array.rows();++i)
        res.push_back(array.row(i));
    return res;
}

bezier_t* wrapBezierConstructor(const point_list_t& array)
{
    t_point_t asVector = vectorFromEigenArray(array);
    return new bezier_t(asVector.begin(), asVector.end());
}


bezier_t* wrapBezierConstructorBounds(const point_list_t& array, const real lb, const real ub)
{
    t_point_t asVector = vectorFromEigenArray(array);
    return new bezier_t(asVector.begin(), asVector.end(), lb, ub);
}

BOOST_PYTHON_MODULE(spline)
{
    /** BEGIN eigenpy init**/
    eigenpy::enableEigenPy();

    eigenpy::enableEigenPySpecific<point_t,point_t>();
    eigenpy::enableEigenPySpecific<point_list_t,point_list_t>();
    /*eigenpy::exposeAngleAxis();
    eigenpy::exposeQuaternion();*/


    class_<bezier_t>
        ("bezier", no_init)
            .def("__init__", make_constructor(&wrapBezierConstructor))
            .def("__init__", make_constructor(&wrapBezierConstructorBounds))
            .def("min", &bezier_t::min)
            .def("max", &bezier_t::max)
            .def("__call__", &bezier_t::operator())
        ;


    /** END eigenpy init**/

    /** BEGIN enum types **/
  /*#ifdef CLP_FOUND
    enum_<SolverLP>("SolverLP")
            .value("SOLVER_LP_QPOASES", SOLVER_LP_QPOASES)
            .value("SOLVER_LP_CLP", SOLVER_LP_CLP)
            .export_values();
  #else
    enum_<SolverLP>("SolverLP")
            .value("SOLVER_LP_QPOASES", SOLVER_LP_QPOASES)
            .export_values();
  #endif


    enum_<EquilibriumAlgorithm>("EquilibriumAlgorithm")
            .value("EQUILIBRIUM_ALGORITHM_LP", EQUILIBRIUM_ALGORITHM_LP)
            .value("EQUILIBRIUM_ALGORITHM_LP2", EQUILIBRIUM_ALGORITHM_LP2)
            .value("EQUILIBRIUM_ALGORITHM_DLP", EQUILIBRIUM_ALGORITHM_DLP)
            .value("EQUILIBRIUM_ALGORITHM_PP", EQUILIBRIUM_ALGORITHM_PP)
            .value("EQUILIBRIUM_ALGORITHM_IP", EQUILIBRIUM_ALGORITHM_IP)
            .value("EQUILIBRIUM_ALGORITHM_DIP", EQUILIBRIUM_ALGORITHM_DIP)
            .export_values();

    enum_<LP_status>("LP_status")
            .value("LP_STATUS_UNKNOWN", LP_STATUS_UNKNOWN)
            .value("LP_STATUS_OPTIMAL", LP_STATUS_OPTIMAL)
            .value("LP_STATUS_INFEASIBLE", LP_STATUS_INFEASIBLE)
            .value("LP_STATUS_UNBOUNDED", LP_STATUS_UNBOUNDED)
            .value("LP_STATUS_MAX_ITER_REACHED", LP_STATUS_MAX_ITER_REACHED)
            .value("LP_STATUS_ERROR", LP_STATUS_ERROR)
            .export_values();*/

    /** END enum types **/

    /*bool (Equilibrium::*setNewContacts)
            (const MatrixX3ColMajor&, const MatrixX3ColMajor&, const double, const EquilibriumAlgorithm) = &Equilibrium::setNewContacts;

    class_<Equilibrium>("Equilibrium", init<std::string, double, unsigned int, optional <SolverLP, bool, const unsigned int, const bool> >())
            .def("useWarmStart", &Equilibrium::useWarmStart)
            .def("setUseWarmStart", &Equilibrium::setUseWarmStart)
            .def("getName", &Equilibrium::getName)
            .def("getAlgorithm", &Equilibrium::getAlgorithm)
            .def("setNewContacts", setNewContacts)
            .def("computeEquilibriumRobustness", wrapComputeQuasiEquilibriumRobustness)
            .def("computeEquilibriumRobustness", wrapComputeEquilibriumRobustness)
            .def("getPolytopeInequalities", wrapGetPolytopeInequalities)
    ;*/
}

} // namespace spline
