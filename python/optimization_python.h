#include "python_variables.h"
#include "namespace.h"

#include "curves/optimization/definitions.h"
#include "curves/optimization/quadratic_problem.h"

#include <boost/python.hpp>

#ifndef _OPTIMIZATION_PYTHON
#define _OPTIMIZATION_PYTHON

namespace curves
{
namespace optimization
{
namespace python
{

void exposeOptimization();
} // namespace python
} // namespace optimization
} // namespace curves

#endif //_OPTIMIZATION_PYTHON
