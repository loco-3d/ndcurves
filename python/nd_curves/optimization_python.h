#include "python_variables.h"
#include "namespace.h"

#include "nd_curves/optimization/definitions.h"
#include "nd_curves/optimization/quadratic_problem.h"

#include <boost/python.hpp>

#ifndef _OPTIMIZATION_PYTHON
#define _OPTIMIZATION_PYTHON

namespace nd_curves {
namespace optimization {
namespace python {

void exposeOptimization();
}  // namespace python
}  // namespace optimization
}  // namespace nd_curves

#endif  //_OPTIMIZATION_PYTHON
