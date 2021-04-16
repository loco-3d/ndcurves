//
// Copyright (c) 2019 INRIA
//

#ifndef __python_namespace_h__
#define __python_namespace_h__

// Silence a warning about a deprecated use of boost bind by boost python
// at least fo boost 1.73 to 1.75
// ref. https://github.com/stack-of-tasks/tsid/issues/128
#define BOOST_BIND_GLOBAL_PLACEHOLDERS
#include <boost/python.hpp>
#undef BOOST_BIND_GLOBAL_PLACEHOLDERS

namespace ndcurves {
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
}  // namespace ndcurves

#endif  // ifndef __python_namespace_h__
