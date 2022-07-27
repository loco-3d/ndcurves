cmake_minimum_required(VERSION 3.1)

# Project properties
set(PROJECT_ORG loco-3d)
set(PROJECT_NAME ndcurves)
set(PROJECT_DESCRIPTION "creatie and manipulate spline and bezier curves.")
set(PROJECT_URL "https://github.com/${PROJECT_ORG}/${PROJECT_NAME}")

# Project options
option(BUILD_PYTHON_INTERFACE "Build the python bindings" ON)
option(INSTALL_PYTHON_INTERFACE_ONLY "Install *ONLY* the python bindings" OFF)
option(SUFFIX_SO_VERSION "Suffix library name with its version" ON)

# Project configuration
if(NOT INSTALL_PYTHON_INTERFACE_ONLY)
  set(PROJECT_USE_CMAKE_EXPORT TRUE)
endif(NOT INSTALL_PYTHON_INTERFACE_ONLY)
set(CXX_DISABLE_WERROR TRUE)

# JRL-cmakemodule setup
include(cmake/base.cmake)
include(cmake/boost.cmake)

# Project definition
compute_project_args(PROJECT_ARGS LANGUAGES CXX)
project(${PROJECT_NAME} ${PROJECT_ARGS})
check_minimal_cxx_standard(11 ENFORCE)

# Project dependencies
find_package(Eigen3 QUIET NO_CMAKE_PACKAGE_REGISTRY)
if(Eigen3_FOUND)
  add_project_dependency(Eigen3 REQUIRED NO_CMAKE_PACKAGE_REGISTRY
                         PKG_CONFIG_REQUIRES eigen3)
else(Eigen3_FOUND)
  add_project_dependency(Eigen3 MODULE REQUIRED PKG_CONFIG_REQUIRES eigen3)
endif(Eigen3_FOUND)
find_package(pinocchio QUIET)
option(CURVES_WITH_PINOCCHIO_SUPPORT "Build with pinocchio support"
       ${pinocchio_FOUND})
if(CURVES_WITH_PINOCCHIO_SUPPORT)
  if(NOT pinocchio_FOUND)
    message(
      FATAL_ERROR
        "CURVES_WITH_PINOCCHIO_SUPPORT selected, but pinocchio has not been found"
    )
  endif(NOT pinocchio_FOUND)
  add_project_dependency(pinocchio REQUIRED)
  pkg_config_append_cflags("-DCURVES_WITH_PINOCCHIO_SUPPORT")
endif(CURVES_WITH_PINOCCHIO_SUPPORT)
set(PACKAGE_EXTRA_MACROS
    "SET(CURVES_WITH_PINOCCHIO_SUPPORT ${CURVES_WITH_PINOCCHIO_SUPPORT})")

add_project_dependency(Boost REQUIRED COMPONENTS serialization)

if(BUILD_PYTHON_INTERFACE)
  add_project_dependency(eigenpy 2.7.10 REQUIRED PKG_CONFIG_REQUIRES
                         "eigenpy >= 2.7.10")
endif(BUILD_PYTHON_INTERFACE)

# Main Library
set(${PROJECT_NAME}_HEADERS
    include/${PROJECT_NAME}/bernstein.h
    include/${PROJECT_NAME}/bezier_curve.h
    include/${PROJECT_NAME}/constant_curve.h
    include/${PROJECT_NAME}/cross_implementation.h
    include/${PROJECT_NAME}/cubic_hermite_spline.h
    include/${PROJECT_NAME}/curve_abc.h
    include/${PROJECT_NAME}/curve_constraint.h
    include/${PROJECT_NAME}/curve_conversion.h
    include/${PROJECT_NAME}/exact_cubic.h
    include/${PROJECT_NAME}/fwd.h
    include/${PROJECT_NAME}/helpers/effector_spline.h
    include/${PROJECT_NAME}/helpers/effector_spline_rotation.h
    include/${PROJECT_NAME}/linear_variable.h
    include/${PROJECT_NAME}/MathDefs.h
    include/${PROJECT_NAME}/optimization/definitions.h
    include/${PROJECT_NAME}/optimization/details.h
    include/${PROJECT_NAME}/optimization/integral_cost.h
    include/${PROJECT_NAME}/optimization/quadratic_problem.h
    include/${PROJECT_NAME}/piecewise_curve.h
    include/${PROJECT_NAME}/polynomial.h
    include/${PROJECT_NAME}/python/python_definitions.h
    include/${PROJECT_NAME}/quadratic_variable.h
    include/${PROJECT_NAME}/se3_curve.h
    include/${PROJECT_NAME}/serialization/archive.hpp
    include/${PROJECT_NAME}/serialization/curves.hpp
    include/${PROJECT_NAME}/serialization/eigen-matrix.hpp
    include/${PROJECT_NAME}/serialization/registeration.hpp
    include/${PROJECT_NAME}/sinusoidal.h
    include/${PROJECT_NAME}/so3_linear.h)

add_library(${PROJECT_NAME} INTERFACE)
modernize_target_link_libraries(
  ${PROJECT_NAME}
  SCOPE
  INTERFACE
  TARGETS
  Eigen3::Eigen
  INCLUDE_DIRS
  ${EIGEN3_INCLUDE_DIR})
target_include_directories(${PROJECT_NAME}
                           INTERFACE $<INSTALL_INTERFACE:include>)
target_link_libraries(${PROJECT_NAME} INTERFACE Boost::serialization)
if(CURVES_WITH_PINOCCHIO_SUPPORT)
  target_link_libraries(${PROJECT_NAME} INTERFACE pinocchio::pinocchio)
  target_compile_definitions(${PROJECT_NAME}
                             INTERFACE -DCURVES_WITH_PINOCCHIO_SUPPORT)
endif(CURVES_WITH_PINOCCHIO_SUPPORT)

if(NOT INSTALL_PYTHON_INTERFACE_ONLY)
  install(
    TARGETS ${PROJECT_NAME}
    EXPORT ${TARGETS_EXPORT_NAME}
    DESTINATION lib)
  install(FILES package.xml DESTINATION share/${PROJECT_NAME})
endif(NOT INSTALL_PYTHON_INTERFACE_ONLY)

if(BUILD_PYTHON_INTERFACE)
  add_subdirectory(python)
endif(BUILD_PYTHON_INTERFACE)

add_subdirectory(tests)
