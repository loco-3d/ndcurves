ADD_SUBDIRECTORY(${PY_NAME})

PYTHON_INSTALL_ON_SITE(${PY_NAME} __init__.py)
PYTHON_INSTALL_ON_SITE(${PY_NAME} plot.py)
PYTHON_INSTALL_ON_SITE(${PY_NAME} optimization.py)

ADD_PYTHON_UNIT_TEST("python-curves" "python/test/test.py" "python")
ADD_PYTHON_UNIT_TEST("python-optimization" "python/test/optimization.py" "python")
ADD_PYTHON_UNIT_TEST("python-notebook" "python/test/notebook.py" "python")
