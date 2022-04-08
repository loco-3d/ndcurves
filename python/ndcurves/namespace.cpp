//
// Copyright (c) 2019 INRIA
//

#include "namespace.h"

namespace ndcurves {
namespace python {

bp::object getOrCreatePythonNamespace(const std::string& submodule_name) {
  bp::scope current_scope;
  std::string current_scope_name(
      bp::extract<const char*>(current_scope.attr("__name__")));
  std::string complete_submodule_name =
      current_scope_name + "." + submodule_name;

  bp::object submodule(
      bp::borrowed(PyImport_AddModule(complete_submodule_name.c_str())));
  current_scope.attr(submodule_name.c_str()) = submodule;

  return submodule;
}
}  // namespace python
}  // namespace ndcurves
