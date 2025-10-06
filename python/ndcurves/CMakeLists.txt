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
set(${wrap}_INSTALL_DIR ${PYTHON_SITELIB}/${PROJECT_NAME})
add_library(${wrap} SHARED ${${PROJECT_NAME}_WRAP_SOURCES})
get_relative_rpath(${${wrap}_INSTALL_DIR} ${wrap}_INSTALL_RPATH)
set_target_properties(
  ${wrap}
  PROPERTIES PREFIX ""
             SUFFIX ${PYTHON_EXT_SUFFIX}
             OUTPUT_NAME ${PROJECT_NAME}
             INSTALL_RPATH "${${wrap}_INSTALL_RPATH}")
cxx_flags_by_compiler_frontend(GNU -Wno-conversion OUTPUT PRIVATE_OPTIONS
                               FILTER)
target_compile_options(${wrap} PRIVATE ${PRIVATE_OPTIONS})
target_link_libraries(${wrap} PUBLIC ${PROJECT_NAME}::${PROJECT_NAME}
                                     eigenpy::eigenpy)

if(GENERATE_PYTHON_STUBS)
  load_stubgen()
  generate_stubs(${CMAKE_CURRENT_BINARY_DIR} ${PROJECT_NAME} ${PYTHON_SITELIB}
                 ${wrap})
endif(GENERATE_PYTHON_STUBS)

if(APPLE)
  # We need to change the extension for python bindings
  set_target_properties(${wrap} PROPERTIES SUFFIX ".so")
endif(APPLE)

install(TARGETS ${wrap} DESTINATION ${${wrap}_INSTALL_DIR})
