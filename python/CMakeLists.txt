cmake_minimum_required( VERSION 2.8 )

include_directories("${EIGEN3_INCLUDE_DIR}")

ADD_REQUIRED_DEPENDENCY("eigenpy")

include_directories("${PROJECT_SOURCE_DIR}/include")
FILE(GLOB_RECURSE HeaderFiles "${PROJECT_SOURCE_DIR}/include/spline/*.h")

# Define the wrapper library that wraps our library
add_library( spline SHARED spline_python.cpp )
#~ target_link_libraries( centroidal_dynamics ${Boost_LIBRARIES} ${PROJECT_NAME} )
# don't prepend wrapper library name with lib
set_target_properties( spline PROPERTIES PREFIX "" )

IF(APPLE)
	# We need to change the extension for python bindings
	SET_TARGET_PROPERTIES(spline PROPERTIES SUFFIX ".so")
ENDIF(APPLE)
PKG_CONFIG_USE_DEPENDENCY(spline eigenpy)

INSTALL(
	FILES ${LIBRARY_OUTPUT_PATH}spline.so	DESTINATION ${PYTHON_SITELIB}
)

ADD_PYTHON_UNIT_TEST("python-spline" "python/test/test.py" "python")
