#include "python_variables.h"
#include "archive_python_binding.h"
#include "namespace.h"
#include "curves/optimization/definitions.h"

#include <boost/python.hpp>
#include <boost/python/enum.hpp>

namespace curves
{
namespace optimization
{
  namespace python
  {
    void exposeOptimization()
    {
        // using the optimization scope
        bp::scope current_scope = curves::python::getOrCreatePythonNamespace("optimization");
        /** BEGIN enums**/
        bp::enum_<constraint_flag>("constraint_flag")
                .value("INIT_POS" , INIT_POS)
                .value("INIT_POS" , INIT_VEL)
                .value("INIT_ACC" , INIT_ACC)
                .value("INIT_JERK", INIT_JERK)
                .value("END_POS"  , END_POS)
                .value("END_POS"  , END_VEL)
                .value("END_ACC"  , END_ACC)
                .value("END_JERK" , END_JERK)
                .export_values();
        /** END curve constraints**/
    }

  } // namespace python
} // namespace optimization
} // namespace curves
