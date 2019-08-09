// Copyright (c) 2015-2018, CNRS
// Authors: Justin Carpentier <jcarpent@laas.fr>

#ifndef __curves_python_serialization_archive_hpp__
#define __curves_python_serialization_archive_hpp__

#include <boost/python.hpp>
#include "curves/curve_abc.h"

namespace curves
{
  namespace bp = boost::python;
  template<typename Derived>
  struct SerializableVisitor
    : public boost::python::def_visitor< SerializableVisitor<Derived> >
  {

    // TO DO !!!!! Try to fix and remove all the .def for serialization in curves_python

    template<class PyClass>
    void visit(PyClass& cl) const
    {
      cl
      /*
      .def("saveAsText",&Derived::saveAsText<Derived>,bp::args("filename"),"Saves *this inside a text file.")
      .def("loadFromText",Derived::loadFromText<Derived>,bp::args("filename"),"Loads *this from a text file.")
      .def("saveAsXML",Derived::saveAsXML<Derived>,bp::args("filename","tag_name"),"Saves *this inside a XML file.")
      .def("loadFromXML",Derived::loadFromXML<Derived>,bp::args("filename","tag_name"),"Loads *this from a XML file.")
      .def("saveAsBinary",Derived::saveAsBinary<Derived>,bp::args("filename"),"Saves *this inside a binary file.")
      .def("loadFromBinary",Derived::loadFromBinary<Derived>,bp::args("filename"),"Loads *this from a binary file.")
      */
      ;
    }

  };
}

#endif // ifndef __multicontact_api_python_serialization_archive_hpp__
