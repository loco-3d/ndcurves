CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

INCLUDE(cmake/base.cmake)
INCLUDE(cmake/boost.cmake)
INCLUDE(cmake/test.cmake)
INCLUDE(cmake/python.cmake)
INCLUDE(cmake/hpp.cmake)

SET(PROJECT_ORG loco-3d)
SET(PROJECT_NAME curves)
SET(PROJECT_DESCRIPTION
  "template based classes for creating and manipulating spline and bezier curves. Comes with extra options specific to end-effector trajectories in robotics."
  )

# Disable -Werror on Unix for now.
SET(CXX_DISABLE_WERROR True)

SETUP_HPP_PROJECT()

ADD_REQUIRED_DEPENDENCY(eigen3)

SET(BOOST_COMPONENTS unit_test_framework)

OPTION (BUILD_PYTHON_INTERFACE "Build the python binding" ON)
IF(BUILD_PYTHON_INTERFACE)
  FINDPYTHON()
  INCLUDE_DIRECTORIES(SYSTEM ${PYTHON_INCLUDE_DIRS})

  STRING(REGEX REPLACE "-" "_" PY_NAME ${PROJECT_NAME})

  ADD_REQUIRED_DEPENDENCY("eigenpy")
  SET(BOOST_COMPONENTS ${BOOST_COMPONENTS} python)

  ADD_SUBDIRECTORY(python)
ENDIF(BUILD_PYTHON_INTERFACE)

SEARCH_FOR_BOOST()
INCLUDE_DIRECTORIES(SYSTEM ${Boost_INCLUDE_DIRS})

ADD_SUBDIRECTORY(include/curve)
ADD_SUBDIRECTORY(tests)

SETUP_HPP_PROJECT_FINALIZE()
