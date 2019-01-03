cmake_minimum_required(VERSION 2.6)
project(spline)
INCLUDE(cmake/base.cmake)
INCLUDE(cmake/test.cmake)
INCLUDE(cmake/python.cmake)

SET(PROJECT_ORG humanoid-path-planner)
SET(PROJECT_NAME hpp-spline)
SET(PROJECT_DESCRIPTION
	 "template based classes for creating and manipulating spline and bezier curves. Comes with extra options specific to end-effector trajectories in robotics."
 )
SET(PROJECT_URL "https://github.com/${PROJECT_ORG}/${PROJECT_NAME}")

# Disable -Werror on Unix for now.
SET(CXX_DISABLE_WERROR True)
SET(CMAKE_VERBOSE_MAKEFILE True)

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

SETUP_PROJECT()

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

add_subdirectory (tests)

install(FILES
	${CMAKE_SOURCE_DIR}/include/spline/bernstein.h
	${CMAKE_SOURCE_DIR}/include/spline/bezier_polynom_conversion.h
	${CMAKE_SOURCE_DIR}/include/spline/curve_abc.h
	${CMAKE_SOURCE_DIR}/include/spline/exact_cubic.h
	${CMAKE_SOURCE_DIR}/include/spline/MathDefs.h
	${CMAKE_SOURCE_DIR}/include/spline/polynom.h
	${CMAKE_SOURCE_DIR}/include/spline/spline_deriv_constraint.h
	${CMAKE_SOURCE_DIR}/include/spline/bezier_curve.h
	${CMAKE_SOURCE_DIR}/include/spline/cubic_spline.h
	${CMAKE_SOURCE_DIR}/include/spline/curve_constraint.h
	${CMAKE_SOURCE_DIR}/include/spline/quintic_spline.h
	DESTINATION ${CMAKE_INSTALL_PREFIX}/include/spline
	)

install(FILES
	${CMAKE_SOURCE_DIR}/include/spline/helpers/effector_spline.h
	${CMAKE_SOURCE_DIR}/include/spline/helpers/effector_spline_rotation.h
	DESTINATION ${CMAKE_INSTALL_PREFIX}/include/spline/helpers
	)

SETUP_PROJECT_FINALIZE()
