ADD_SUBDIRECTORY(${PROJECT_NAME})

PYTHON_INSTALL_ON_SITE(${PROJECT_NAME} __init__.py)
PYTHON_INSTALL_ON_SITE(${PROJECT_NAME} plot.py)
PYTHON_INSTALL_ON_SITE(${PROJECT_NAME} optimization.py)

ADD_PYTHON_UNIT_TEST("python-curves" "python/test/test.py" "python")
ADD_PYTHON_UNIT_TEST("python-constant" "python/test/test-constant.py" "python")
ADD_PYTHON_UNIT_TEST("python-optimization" "python/test/optimization.py" "python")
ADD_PYTHON_UNIT_TEST("python-notebook" "python/test/notebook.py" "python")
