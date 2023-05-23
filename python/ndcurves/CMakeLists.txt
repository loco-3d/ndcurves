if(GENERATE_PYTHON_STUBS)
  include("${JRL_CMAKE_MODULES}/stubs.cmake")
endif(GENERATE_PYTHON_STUBS)

set(${PROJECT_NAME}_WRAP_SOURCES
    curves_python.cpp
    optimization_python.cpp
    optimization_python.h
    python_variables.cpp
    python_variables.h
    namespace.cpp
    namespace.h)

set(wrap ${PROJECT_NAME}_wrap)
add_library(${wrap} SHARED ${${PROJECT_NAME}_WRAP_SOURCES})
set_target_properties(
  ${wrap}
  PROPERTIES PREFIX ""
             OUTPUT_NAME ${PROJECT_NAME}
             INSTALL_RPATH "\$ORIGIN/../../..")
target_compile_options(${wrap} PRIVATE "-Wno-conversion")
target_link_libraries(${wrap} PUBLIC ${PROJECT_NAME})
target_link_libraries(${wrap} PUBLIC eigenpy::eigenpy)
target_link_boost_python(${wrap} PUBLIC)

if(GENERATE_PYTHON_STUBS)
  load_stubgen()
  generate_stubs(${CMAKE_CURRENT_BINARY_DIR} ${PROJECT_NAME} ${PYTHON_SITELIB}
                 ${wrap})
endif(GENERATE_PYTHON_STUBS)

if(APPLE)
  # We need to change the extension for python bindings
  set_target_properties(${wrap} PROPERTIES SUFFIX ".so")
endif(APPLE)

install(TARGETS ${wrap} DESTINATION ${PYTHON_SITELIB}/${PROJECT_NAME})
