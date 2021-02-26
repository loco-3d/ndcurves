#include "python_variables.h"
#include "namespace.h"

#include "NDcurves/optimization/definitions.h"
#include "NDcurves/optimization/quadratic_problem.h"

#include <boost/python.hpp>

#ifndef _OPTIMIZATION_PYTHON
#define _OPTIMIZATION_PYTHON

namespace NDcurves {
namespace optimization {
namespace python {

void exposeOptimization();
}  // namespace python
}  // namespace optimization
}  // namespace NDcurves

#endif  //_OPTIMIZATION_PYTHON
