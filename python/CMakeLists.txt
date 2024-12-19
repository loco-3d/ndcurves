add_subdirectory(${PROJECT_NAME})

python_install_on_site(${PROJECT_NAME} __init__.py)
python_install_on_site(${PROJECT_NAME} plot.py)
python_install_on_site(${PROJECT_NAME} optimization.py)

if(BUILD_TESTING)
  add_python_unit_test("python-curves" "python/test/test.py" "python")
  add_python_unit_test("python-constant" "python/test/test-constant.py"
                       "python")
  add_python_unit_test("python-sinusoidal" "python/test/test-sinusoidal.py"
                       "python")
  add_python_unit_test("python-minjerk" "python/test/test-minjerk.py" "python")
  add_python_unit_test("python-optimization" "python/test/optimization.py"
                       "python")
  add_python_unit_test("python-notebook" "python/test/notebook.py" "python")
  add_python_unit_test("python-registration" "python/test/registration.py"
                       "python")
  add_python_unit_test("python-curve-constraints"
                       "python/test/test-curve-constraints.py" "python")
endif()
