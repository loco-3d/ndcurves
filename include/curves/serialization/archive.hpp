// Copyright (c) 2015-2018, CNRS
// Authors: Justin Carpentier <jcarpent@laas.fr>

#ifndef __curves_serialization_archive_hpp__
#define __curves_serialization_archive_hpp__
#include <fstream>
#include <string>
#include <stdexcept>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include "registeration.hpp"

namespace curves {
namespace serialization {
struct Serializable {
 private:
  template <class Derived>
  Derived& derived() {
    return *static_cast<Derived*>(this);
  }
  template <class Derived>
  const Derived& derived() const {
    return *static_cast<const Derived*>(this);
  }

 public:
  /// \brief Loads a Derived object from a text file.
  template <class Derived>
  void loadFromText(const std::string& filename) {
    std::ifstream ifs(filename.c_str());
    if (ifs) {
      boost::archive::text_iarchive ia(ifs);
      register_types<boost::archive::text_iarchive>(ia);
      ia >> derived<Derived>();
    } else {
      const std::string exception_message(filename + " does not seem to be a valid file.");
      throw std::invalid_argument(exception_message);
    }
  }

  /// \brief Saved a Derived object as a text file.
  template <class Derived>
  void saveAsText(const std::string& filename) const {
    std::ofstream ofs(filename.c_str());
    if (ofs) {
      boost::archive::text_oarchive oa(ofs);
      register_types<boost::archive::text_oarchive>(oa);
      oa << derived<Derived>();
    } else {
      const std::string exception_message(filename + " does not seem to be a valid file.");
      throw std::invalid_argument(exception_message);
    }
  }

  /// \brief Loads a Derived object from an XML file.
  template <class Derived>
  void loadFromXML(const std::string& filename, const std::string& tag_name) {
    if (tag_name.empty()) {
      throw std::invalid_argument("tag_name cannot be empty.");
    }
    std::ifstream ifs(filename.c_str());
    if (ifs) {
      boost::archive::xml_iarchive ia(ifs);
      register_types<boost::archive::xml_iarchive>(ia);
      ia >> boost::serialization::make_nvp(tag_name.c_str(), derived<Derived>());
    } else {
      const std::string exception_message(filename + " does not seem to be a valid file.");
      throw std::invalid_argument(exception_message);
    }
  }

  /// \brief Saved a Derived object as an XML file.
  template <class Derived>
  void saveAsXML(const std::string& filename, const std::string& tag_name) const {
    if (tag_name.empty()) {
      throw std::invalid_argument("tag_name cannot be empty.");
    }
    std::ofstream ofs(filename.c_str());
    if (ofs) {
      boost::archive::xml_oarchive oa(ofs);
      register_types<boost::archive::xml_oarchive>(oa);
      oa << boost::serialization::make_nvp(tag_name.c_str(), derived<Derived>());
    } else {
      const std::string exception_message(filename + " does not seem to be a valid file.");
      throw std::invalid_argument(exception_message);
    }
  }

  /// \brief Loads a Derived object from an binary file.
  template <class Derived>
  void loadFromBinary(const std::string& filename) {
    std::ifstream ifs(filename.c_str());
    if (ifs) {
      boost::archive::binary_iarchive ia(ifs);
      register_types<boost::archive::binary_iarchive>(ia);
      ia >> derived<Derived>();
    } else {
      const std::string exception_message(filename + " does not seem to be a valid file.");
      throw std::invalid_argument(exception_message);
    }
  }

  /// \brief Saved a Derived object as an binary file.
  template <class Derived>
  void saveAsBinary(const std::string& filename) const {
    std::ofstream ofs(filename.c_str());
    if (ofs) {
      boost::archive::binary_oarchive oa(ofs);
      register_types<boost::archive::binary_oarchive>(oa);
      oa << derived<Derived>();
    } else {
      const std::string exception_message(filename + " does not seem to be a valid file.");
      throw std::invalid_argument(exception_message);
    }
  }
};  // End struct Serializable

}  // namespace serialization

}  // namespace curves

#endif  // ifndef __multicontact_api_serialization_archive_hpp__
