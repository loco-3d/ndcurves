CMAKE_MINIMUM_REQUIRED(VERSION 3.1)

# Project properties
SET(PROJECT_ORG loco-3d)
SET(PROJECT_NAME curves)
SET(PROJECT_DESCRIPTION "creatie and manipulate spline and bezier curves.")
SET(PROJECT_URL "https://github.com/${PROJECT_ORG}/${PROJECT_NAME}")

# Project options
OPTION(BUILD_PYTHON_INTERFACE "Build the python bindings" ON)
OPTION(INSTALL_PYTHON_INTERFACE_ONLY "Install *ONLY* the python bindings" OFF)
OPTION(SUFFIX_SO_VERSION "Suffix library name with its version" ON)

# Project configuration
IF(NOT INSTALL_PYTHON_INTERFACE_ONLY)
  SET(PROJECT_USE_CMAKE_EXPORT TRUE)
ENDIF(NOT INSTALL_PYTHON_INTERFACE_ONLY)
SET(CXX_DISABLE_WERROR TRUE)

# JRL-cmakemodule setup
INCLUDE(cmake/base.cmake)
INCLUDE(cmake/boost.cmake)
INCLUDE(cmake/python.cmake)

# Project definition
COMPUTE_PROJECT_ARGS(PROJECT_ARGS LANGUAGES CXX)
PROJECT(${PROJECT_NAME} ${PROJECT_ARGS})

# Project dependencies
ADD_PROJECT_DEPENDENCY(Eigen3 REQUIRED PKG_CONFIG_REQUIRES eigen3)
ADD_PROJECT_DEPENDENCY(pinocchio)
OPTION(CURVES_WITH_PINOCCHIO_SUPPORT "Build with pinocchio support" ${pinocchio_FOUND})
IF(CURVES_WITH_PINOCCHIO_SUPPORT)
  IF(NOT pinocchio_FOUND)
    MESSAGE(FATAL_ERROR "CURVES_WITH_PINOCCHIO_SUPPORT selected, but pinocchio has not been found")
  ENDIF(NOT pinocchio_FOUND)
  ADD_DEFINITIONS(-DCURVES_WITH_PINOCCHIO_SUPPORT)
  PKG_CONFIG_APPEND_CFLAGS("-DCURVES_WITH_PINOCCHIO_SUPPORT")
ENDIF(CURVES_WITH_PINOCCHIO_SUPPORT)
SET(PACKAGE_EXTRA_MACROS "SET(CURVES_WITH_PINOCCHIO_SUPPORT ${CURVES_WITH_PINOCCHIO_SUPPORT})")

SET(BOOST_COMPONENTS unit_test_framework serialization)

IF(BUILD_PYTHON_INTERFACE)
  FINDPYTHON()
  ADD_PROJECT_DEPENDENCY(eigenpy REQUIRED PKG_CONFIG_REQUIRES "eigenpy >= 1.6.8")
  SET(BOOST_COMPONENTS ${BOOST_COMPONENTS} python)
ENDIF(BUILD_PYTHON_INTERFACE)

SEARCH_FOR_BOOST()

# Main Library
SET(${PROJECT_NAME}_HEADERS
  include/${PROJECT_NAME}/bernstein.h
  include/${PROJECT_NAME}/curve_conversion.h
  include/${PROJECT_NAME}/curve_abc.h
  include/${PROJECT_NAME}/exact_cubic.h
  include/${PROJECT_NAME}/MathDefs.h
  include/${PROJECT_NAME}/polynomial.h
  include/${PROJECT_NAME}/bezier_curve.h
  include/${PROJECT_NAME}/cubic_spline.h
  include/${PROJECT_NAME}/curve_constraint.h
  include/${PROJECT_NAME}/quintic_spline.h
  include/${PROJECT_NAME}/linear_variable.h
  include/${PROJECT_NAME}/quadratic_variable.h
  include/${PROJECT_NAME}/cubic_hermite_spline.h
  include/${PROJECT_NAME}/piecewise_curve.h
  include/${PROJECT_NAME}/so3_linear.h
  include/${PROJECT_NAME}/se3_curve.h
  include/${PROJECT_NAME}/sinusoidal.h
  include/${PROJECT_NAME}/fwd.h
  include/${PROJECT_NAME}/constant_curve.h
  include/${PROJECT_NAME}/helpers/effector_spline.h
  include/${PROJECT_NAME}/helpers/effector_spline_rotation.h
  include/${PROJECT_NAME}/optimization/definitions.h
  include/${PROJECT_NAME}/optimization/details.h
  include/${PROJECT_NAME}/optimization/quadratic_problem.h
  include/${PROJECT_NAME}/optimization/integral_cost.h
  include/${PROJECT_NAME}/python/python_definitions.h
  include/${PROJECT_NAME}/serialization/archive.hpp
  include/${PROJECT_NAME}/serialization/registeration.hpp
  include/${PROJECT_NAME}/serialization/eigen-matrix.hpp
  include/${PROJECT_NAME}/serialization/curves.hpp
  )

ADD_LIBRARY(${PROJECT_NAME} INTERFACE)
TARGET_INCLUDE_DIRECTORIES(${PROJECT_NAME} SYSTEM INTERFACE ${EIGEN3_INCLUDE_DIRS})
TARGET_INCLUDE_DIRECTORIES(${PROJECT_NAME} INTERFACE $<INSTALL_INTERFACE:include>)
IF(CURVES_WITH_PINOCCHIO_SUPPORT)
  TARGET_LINK_LIBRARIES(${PROJECT_NAME} INTERFACE pinocchio::pinocchio)
ENDIF(CURVES_WITH_PINOCCHIO_SUPPORT)

IF(NOT INSTALL_PYTHON_INTERFACE_ONLY)
  INSTALL(TARGETS ${PROJECT_NAME} EXPORT ${TARGETS_EXPORT_NAME} DESTINATION lib)
ENDIF(NOT INSTALL_PYTHON_INTERFACE_ONLY)

IF(BUILD_PYTHON_INTERFACE)
  ADD_SUBDIRECTORY(python)
ENDIF(BUILD_PYTHON_INTERFACE)

ADD_SUBDIRECTORY(tests)
