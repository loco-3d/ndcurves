cmake_minimum_required(VERSION 3.22)

# Project properties
set(PROJECT_ORG loco-3d)
set(PROJECT_NAME ndcurves)
set(PROJECT_DESCRIPTION "create and manipulate spline and bezier curves.")
set(PROJECT_URL "https://github.com/${PROJECT_ORG}/${PROJECT_NAME}")

# Project options
option(SUFFIX_SO_VERSION "Suffix library name with its version" ON)

# Project configuration
set(PROJECT_USE_CMAKE_EXPORT TRUE)
set(CXX_DISABLE_WERROR TRUE)

# JRL-cmakemodule setup
set(JRL_CMAKE_MODULES "${CMAKE_CURRENT_LIST_DIR}/cmake")
if(EXISTS "${JRL_CMAKE_MODULES}/base.cmake")
  message(STATUS "JRL cmakemodules found in 'cmake/' git submodule")
else()
  find_package(jrl-cmakemodules QUIET CONFIG)
  if(jrl-cmakemodules_FOUND)
    get_property(
      JRL_CMAKE_MODULES
      TARGET jrl-cmakemodules::jrl-cmakemodules
      PROPERTY INTERFACE_INCLUDE_DIRECTORIES)
    message(STATUS "JRL cmakemodules found on system at ${JRL_CMAKE_MODULES}")
  else()
    message(STATUS "JRL cmakemodules not found. Let's fetch it.")
    include(FetchContent)
    FetchContent_Declare(
      "jrl-cmakemodules"
      GIT_REPOSITORY "https://github.com/jrl-umi3218/jrl-cmakemodules.git")
    FetchContent_MakeAvailable("jrl-cmakemodules")
    FetchContent_GetProperties("jrl-cmakemodules" SOURCE_DIR JRL_CMAKE_MODULES)
  endif()
endif()

include("${JRL_CMAKE_MODULES}/base.cmake")
include("${JRL_CMAKE_MODULES}/boost.cmake")
include(CMakeDependentOption)

cmake_dependent_option(
  GENERATE_PYTHON_STUBS
  "Generate the Python stubs associated to the Python library" OFF
  BUILD_PYTHON_INTERFACE ON)

set_default_cmake_build_type("Release")

# Project definition
compute_project_args(PROJECT_ARGS LANGUAGES CXX)
project(${PROJECT_NAME} ${PROJECT_ARGS})
check_minimal_cxx_standard(14 ENFORCE)

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
  add_project_dependency(pinocchio REQUIRED PKG_CONFIG_REQUIRES)
  pkg_config_append_cflags("-DCURVES_WITH_PINOCCHIO_SUPPORT")
endif(CURVES_WITH_PINOCCHIO_SUPPORT)
set(PACKAGE_EXTRA_MACROS
    "SET(CURVES_WITH_PINOCCHIO_SUPPORT ${CURVES_WITH_PINOCCHIO_SUPPORT})")

add_project_dependency(Boost REQUIRED COMPONENTS serialization)

if(BUILD_PYTHON_INTERFACE)
  add_project_dependency(eigenpy 3.0.0 REQUIRED PKG_CONFIG_REQUIRES
                         "eigenpy >= 3.0.0")
endif()

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
    include/${PROJECT_NAME}/so3_linear.h
    include/${PROJECT_NAME}/so3_smooth.h)

if(BUILD_STANDALONE_PYTHON_INTERFACE)
  add_project_dependency(${PROJECT_NAME} REQUIRED)
  set(PROJECT_EXPORT_NO_TARGET ON)
else()
  add_library(${PROJECT_NAME} INTERFACE)
  add_library(${PROJECT_NAME}::${PROJECT_NAME} ALIAS ${PROJECT_NAME})
  modernize_target_link_libraries(
    ${PROJECT_NAME}
    SCOPE
    INTERFACE
    TARGETS
    Eigen3::Eigen
    INCLUDE_DIRS
    ${EIGEN3_INCLUDE_DIR})
  target_include_directories(
    ${PROJECT_NAME} INTERFACE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
                              $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
  target_link_libraries(${PROJECT_NAME} INTERFACE Boost::serialization)
  if(CURVES_WITH_PINOCCHIO_SUPPORT)
    target_link_libraries(${PROJECT_NAME} INTERFACE pinocchio::pinocchio)
    target_compile_definitions(${PROJECT_NAME}
                               INTERFACE -DCURVES_WITH_PINOCCHIO_SUPPORT)
  endif(CURVES_WITH_PINOCCHIO_SUPPORT)

  # Define M_PI constant in <cmath>
  cxx_flags_by_compiler_frontend(MSVC _USE_MATH_DEFINES OUTPUT
                                 PUBLIC_DEFINITION)
  target_compile_definitions(${PROJECT_NAME} INTERFACE ${PUBLIC_DEFINITION})
endif()

if(NOT BUILD_STANDALONE_PYTHON_INTERFACE)
  install(
    TARGETS ${PROJECT_NAME}
    EXPORT ${TARGETS_EXPORT_NAME}
    DESTINATION lib)
endif()

if(BUILD_PYTHON_INTERFACE)
  add_subdirectory(python)
endif(BUILD_PYTHON_INTERFACE)

if(BUILD_TESTING AND NOT BUILD_STANDALONE_PYTHON_INTERFACE)
  add_subdirectory(tests)
endif()
