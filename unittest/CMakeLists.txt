if(BUILD_PYTHON_INTERFACE)
  add_subdirectory(bindings)
endif(BUILD_PYTHON_INTERFACE)

set(${PROJECT_NAME}_TESTS
    test_curves
    test_constant
    test_sinusoidal
    test_serialization
    test_minjerk
    test_operations
    test_constraints)

macro(ADD_TEST_CFLAGS target flag)
  set_property(
    TARGET ${target}
    APPEND_STRING
    PROPERTY COMPILE_FLAGS " ${flag}")
endmacro(ADD_TEST_CFLAGS)

foreach(TEST ${${PROJECT_NAME}_TESTS})
  add_unit_test(${TEST} ${TEST}.cpp)

  set_tests_properties(${TEST} PROPERTIES TIMEOUT 7200)
  set_target_properties(${TEST} PROPERTIES LINKER_LANGUAGE CXX)
  add_test_cflags(${TEST} "-DBOOST_TEST_DYN_LINK")

  set(MODULE_NAME "${NAME}Test")
  string(REGEX REPLACE "-" "_" MODULE_NAME ${MODULE_NAME})
  add_test_cflags(${TEST} "-DTEST_MODULE_NAME=${MODULE_NAME}")
  add_test_cflags(${TEST} "-DBOOST_TEST_MODULE=${MODULE_NAME}")

  target_link_libraries(${TEST} ${PROJECT_NAME})
  target_compile_definitions(${TEST} PRIVATE -DTEST_DATA_PATH="${CMAKE_CURRENT_SOURCE_DIR}/data/")
endforeach(TEST ${${PROJECT_NAME}_TESTS})
