set(${PROJECT_NAME}_PYTHON_BINDINGS_TESTS
    curves
    constant
    sinusoidal
    minjerk
    optimization
    notebook
    registration
    constraints)

foreach(TEST ${${PROJECT_NAME}_PYTHON_BINDINGS_TESTS})
  python_build(. "test_${TEST}.py")
  add_python_unit_test("test_pybinds_${TEST}"
                       "unittest/bindings/test_${TEST}.py" python)
  add_custom_target(
    "test_pybinds_${TEST}"
    ${CMAKE_COMMAND} -E env
    PYTHONPATH=${PROJECT_BINARY_DIR}/python/ndcurves:$ENV{PYTHONPATH}
    ${PYTHON_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/test_${TEST}.py")
endforeach(TEST ${${PROJECT_NAME}_PYTHON_BINDINGS_TESTS})
