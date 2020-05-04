//
// Copyright (c) 2019 INRIA
//

#ifndef __curves_python_registration_hpp__
#define __curves_python_registration_hpp__

#include <boost/python.hpp>
#include <boost/python/scope.hpp>

#include <eigenpy/registration.hpp>

namespace curves
{
template<typename T>
inline bool register_symbolic_link_to_registered_type()
{
  namespace bp = boost::python;
  if(eigenpy::check_registration<T>())
  {
    const bp::type_info info = bp::type_id<T>();
    const bp::converter::registration* reg = bp::converter::registry::query(info);
    bp::handle<> class_obj(reg->get_class_object());
    bp::scope().attr(reg->get_class_object()->tp_name) = bp::object(class_obj);
    return true;
  }

  return false;
}
} // namespace curves

#endif // ifndef __curves_python_registration_hpp__

