SET(${PY_NAME}_BINDINGS_SOURCES
  curves_python.cpp
  optimization_python.cpp
  optimization_python.h
  python_variables.cpp
  python_variables.h
  namespace.cpp
  namespace.h
  )

ADD_LIBRARY(${PY_NAME} SHARED ${${PY_NAME}_BINDINGS_SOURCES})
SET_TARGET_PROPERTIES(${PY_NAME} PROPERTIES PREFIX "")
TARGET_COMPILE_OPTIONS(${PY_NAME} PRIVATE "-Wno-conversion")
PKG_CONFIG_USE_DEPENDENCY(${PY_NAME} eigenpy)

TARGET_LINK_LIBRARIES(${PY_NAME} ${Boost_LIBRARIES})
# MESSAGE(WARNING "Boost libraries are : ${Boost_LIBRARIES}")

IF(APPLE)
  # We need to change the extension for python bindings
  SET_TARGET_PROPERTIES(${PY_NAME} PROPERTIES SUFFIX ".so")
ENDIF(APPLE)

INSTALL(TARGETS ${PY_NAME} DESTINATION ${PYTHON_SITELIB})

ADD_PYTHON_UNIT_TEST("python-curves" "python/test/test.py" "python")
ADD_PYTHON_UNIT_TEST("python-optimization" "python/test/optimization.py" "python")
