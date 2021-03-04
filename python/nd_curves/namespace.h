//
// Copyright (c) 2019 INRIA
//

#ifndef __python_namespace_h__
#define __python_namespace_h__

#include <boost/python.hpp>

namespace nd_curves {
namespace python {
namespace bp = boost::python;

///
/// \brief Helper to create or simply return an existing namespace in Python
///
/// \param[in] submodule_name name of the submodule
///
/// \returns The submodule related to the namespace name.
///
bp::object getOrCreatePythonNamespace(const std::string& submodule_name);
}  // namespace python
}  // namespace nd_curves

#endif  // ifndef __python_namespace_h__
