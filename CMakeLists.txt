cmake_minimum_required(VERSION 2.6)
project(curve)
INCLUDE(cmake/base.cmake)
INCLUDE(cmake/test.cmake)
INCLUDE(cmake/python.cmake)
INCLUDE(cmake/hpp.cmake)

SET(PROJECT_NAME curves)
SET(PROJECT_DESCRIPTION
  "template based classes for creating and manipulating spline and bezier curves. Comes with extra options specific to end-effector trajectories in robotics."
  )

# Disable -Werror on Unix for now.
SET(CXX_DISABLE_WERROR True)
SET(CMAKE_VERBOSE_MAKEFILE True)

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

SETUP_HPP_PROJECT()

OPTION (BUILD_PYTHON_INTERFACE "Build the python binding" ON)
IF(BUILD_PYTHON_INTERFACE)
  # search for python
  FINDPYTHON(2.7 REQUIRED)
  find_package( PythonLibs 2.7 REQUIRED )
  include_directories( ${PYTHON_INCLUDE_DIRS} )

  find_package( Boost COMPONENTS python REQUIRED )
  include_directories( ${Boost_INCLUDE_DIR} )

  add_subdirectory (python)
ENDIF(BUILD_PYTHON_INTERFACE)


ADD_SUBDIRECTORY(include/hpp/curve)
ADD_SUBDIRECTORY(tests)

SETUP_HPP_PROJECT_FINALIZE()
