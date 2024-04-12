add_unit_test(curves_tests Main.cpp)
target_link_libraries(curves_tests ${PROJECT_NAME})
target_compile_definitions(
  curves_tests PRIVATE -DTEST_DATA_PATH="${CMAKE_CURRENT_SOURCE_DIR}/data/")

set(${PROJECT_NAME}_TESTS
    test-constant
    test-sinusoidal
    test-sinusoidal-serialization
    test-minjerk
    test-operations
    test-curve-constraints
    test-se3
    test-so3-smooth
    test-polynomial)

foreach(TEST ${${PROJECT_NAME}_TESTS})
  add_unit_test(${TEST} ${TEST}.cpp)
  target_link_libraries(${TEST} ${PROJECT_NAME})
endforeach(TEST ${${PROJECT_NAME}_TESTS})
