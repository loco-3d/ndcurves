#include "python_variables.h"
#include "namespace.h"

#include "ndcurves/optimization/definitions.h"
#include "ndcurves/optimization/quadratic_problem.h"

#include <boost/python.hpp>

#ifndef _OPTIMIZATION_PYTHON
#define _OPTIMIZATION_PYTHON

namespace ndcurves {
namespace optimization {
namespace python {

void exposeOptimization();
}  // namespace python
}  // namespace optimization
}  // namespace ndcurves

#endif  //_OPTIMIZATION_PYTHON
