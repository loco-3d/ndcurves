SET(${PY_NAME}_BINDINGS_SOURCES
  curves_python.cpp
  python_variables.cpp
  python_variables.h
  )

ADD_LIBRARY(${PY_NAME} SHARED ${${PY_NAME}_BINDINGS_SOURCES})
SET_TARGET_PROPERTIES(${PY_NAME} PROPERTIES PREFIX "")
PKG_CONFIG_USE_DEPENDENCY(${PY_NAME} eigenpy)

TARGET_LINK_LIBRARIES(${PY_NAME} ${Boost_LIBRARIES})

IF(APPLE)
  # We need to change the extension for python bindings
  SET_TARGET_PROPERTIES(${PY_NAME} PROPERTIES SUFFIX ".so")
ENDIF(APPLE)

INSTALL(TARGETS ${PY_NAME} DESTINATION ${PYTHON_SITELIB})

ADD_PYTHON_UNIT_TEST("python-curves" "python/test/test.py" "python")
