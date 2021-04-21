SET(${PROJECT_NAME}_WRAP_SOURCES
  curves_python.cpp
  optimization_python.cpp
  optimization_python.h
  python_variables.cpp
  python_variables.h
  namespace.cpp
  namespace.h
  )

SET(wrap ${PROJECT_NAME}_wrap)
ADD_LIBRARY(${wrap} SHARED ${${PROJECT_NAME}_WRAP_SOURCES})
SET_TARGET_PROPERTIES(${wrap} PROPERTIES OUTPUT_NAME ${PROJECT_NAME} PREFIX "")
TARGET_COMPILE_OPTIONS(${wrap} PRIVATE "-Wno-conversion")
TARGET_LINK_LIBRARIES(${wrap} ${PROJECT_NAME})
TARGET_LINK_LIBRARIES(${wrap} eigenpy::eigenpy)
TARGET_LINK_BOOST_PYTHON(${wrap})

IF(APPLE)
  # We need to change the extension for python bindings
  SET_TARGET_PROPERTIES(${wrap} PROPERTIES SUFFIX ".so")
ENDIF(APPLE)

INSTALL(TARGETS ${wrap} DESTINATION ${PYTHON_SITELIB}/${PROJECT_NAME})
