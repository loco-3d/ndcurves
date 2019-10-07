#include "python_variables.h"
#include "namespace.h"

#include <boost/python.hpp>

#ifndef _OPTIMIZATION_PYTHON
#define _VARIABLES_PYTHON_BINDINGS

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
