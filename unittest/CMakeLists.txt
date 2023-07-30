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

foreach(TEST ${${PROJECT_NAME}_TESTS})
  add_unit_test(${TEST} ${TEST}.cpp)
  target_link_libraries(${TEST} ${PROJECT_NAME})
  target_compile_definitions(${TEST} PRIVATE -DTEST_DATA_PATH="${CMAKE_CURRENT_SOURCE_DIR}/data/")
endforeach(TEST ${${PROJECT_NAME}_TESTS})
