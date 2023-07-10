// Copyright (c) 2015-2018, CNRS
// Authors: Justin Carpentier <jcarpent@laas.fr>

#ifndef __curves_python_serialization_archive_hpp__
#define __curves_python_serialization_archive_hpp__

#include <boost/python.hpp>

namespace ndcurves {
namespace bp = boost::python;
template <typename Derived>
struct SerializableVisitor
    : public boost::python::def_visitor<SerializableVisitor<Derived> > {
  template <class PyClass>
  void visit(PyClass& cl) const {
    cl.def("saveAsText", &Derived::template saveAsText<Derived>,
           bp::args("self", "filename"), "Saves *this inside a text file.")
        .def("loadFromText", &Derived::template loadFromText<Derived>,
             bp::args("self", "filename"), "Loads *this from a text file.")
        .def("saveAsXML", &Derived::template saveAsXML<Derived>,
             bp::args("self", "filename", "tag_name"),
             "Saves *this inside a XML file.")
        .def("loadFromXML", &Derived::template loadFromXML<Derived>,
             bp::args("self", "filename", "tag_name"),
             "Loads *this from a XML file.")
        .def("saveAsBinary", &Derived::template saveAsBinary<Derived>,
             bp::args("self", "filename"), "Saves *this inside a binary file.")
        .def("loadFromBinary", &Derived::template loadFromBinary<Derived>,
             bp::args("self", "filename"), "Loads *this from a binary file.");
  }
};
}  // namespace ndcurves

#endif  // ifndef __multicontact_api_python_serialization_archive_hpp__
