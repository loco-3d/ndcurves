ADD_UNIT_TEST(curves_tests Main.cpp)
TARGET_LINK_LIBRARIES(curves_tests ${PROJECT_NAME})
TARGET_COMPILE_DEFINITIONS(curves_tests PRIVATE -DTEST_DATA_PATH="${CMAKE_CURRENT_SOURCE_DIR}/data/")

SET(${PROJECT_NAME}_TESTS
  test-constant
  test-sinusoidal
  test-sinusoidal-serialization
  test-minjerk
  test-operations
  )

FOREACH(TEST ${${PROJECT_NAME}_TESTS})
  ADD_UNIT_TEST(${TEST} ${TEST})
  TARGET_LINK_LIBRARIES(${TEST} ${PROJECT_NAME})
ENDFOREACH(TEST ${${PROJECT_NAME}_TESTS})
