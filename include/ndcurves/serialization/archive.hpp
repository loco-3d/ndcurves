// Copyright (c) 2015-2018, CNRS
// Authors: Justin Carpentier <jcarpent@laas.fr>

#ifndef __curves_serialization_archive_hpp__
#define __curves_serialization_archive_hpp__
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/serialization/version.hpp>
#include <fstream>
#include <stdexcept>
#include <string>

/* Define the current version number for the serialization
 * Must be increased everytime the save() method of a class is modified
 * Or when a change is made to register_types()
 * */
const unsigned int CURVES_API_VERSION = 1;

#define SINGLE_ARG(...) \
  __VA_ARGS__  // Macro used to be able to put comma in the following macro
               // arguments
// Macro used to define the serialization version of a templated class
#define DEFINE_CLASS_TEMPLATE_VERSION(Template, Type)         \
  namespace boost {                                           \
  namespace serialization {                                   \
  template <Template>                                         \
  struct version<Type> {                                      \
    static constexpr unsigned int value = CURVES_API_VERSION; \
  };                                                          \
  template <Template>                                         \
  constexpr unsigned int version<Type>::value;                \
  }                                                           \
  }

namespace ndcurves {
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
      ia >> derived<Derived>();
    } else {
      const std::string exception_message(filename +
                                          " does not seem to be a valid file.");
      throw std::invalid_argument(exception_message);
    }
  }

  /// \brief Saved a Derived object as a text file.
  template <class Derived>
  void saveAsText(const std::string& filename) const {
    std::ofstream ofs(filename.c_str());
    if (ofs) {
      boost::archive::text_oarchive oa(ofs);
      oa << derived<Derived>();
    } else {
      const std::string exception_message(filename +
                                          " does not seem to be a valid file.");
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
      ia >>
          boost::serialization::make_nvp(tag_name.c_str(), derived<Derived>());
    } else {
      const std::string exception_message(filename +
                                          " does not seem to be a valid file.");
      throw std::invalid_argument(exception_message);
    }
  }

  /// \brief Saved a Derived object as an XML file.
  template <class Derived>
  void saveAsXML(const std::string& filename,
                 const std::string& tag_name) const {
    if (tag_name.empty()) {
      throw std::invalid_argument("tag_name cannot be empty.");
    }
    std::ofstream ofs(filename.c_str());
    if (ofs) {
      boost::archive::xml_oarchive oa(ofs);
      oa << boost::serialization::make_nvp(tag_name.c_str(),
                                           derived<Derived>());
    } else {
      const std::string exception_message(filename +
                                          " does not seem to be a valid file.");
      throw std::invalid_argument(exception_message);
    }
  }

  /// \brief Loads a Derived object from an binary file.
  template <class Derived>
  void loadFromBinary(const std::string& filename) {
    std::ifstream ifs(filename.c_str(), std::ios::binary);
    if (ifs) {
      boost::archive::binary_iarchive ia(ifs);
      ia >> derived<Derived>();
    } else {
      const std::string exception_message(filename +
                                          " does not seem to be a valid file.");
      throw std::invalid_argument(exception_message);
    }
  }

  /// \brief Saved a Derived object as an binary file.
  template <class Derived>
  void saveAsBinary(const std::string& filename) const {
    std::ofstream ofs(filename.c_str(), std::ios::binary);
    if (ofs) {
      boost::archive::binary_oarchive oa(ofs);
      oa << derived<Derived>();
    } else {
      const std::string exception_message(filename +
                                          " does not seem to be a valid file.");
      throw std::invalid_argument(exception_message);
    }
  }
};  // End struct Serializable

}  // namespace serialization

}  // namespace ndcurves

#endif  // ifndef __multicontact_api_serialization_archive_hpp__
