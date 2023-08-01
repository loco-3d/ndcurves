add_subdirectory(${PROJECT_NAME})

python_install_on_site(${PROJECT_NAME} __init__.py)
python_install_on_site(${PROJECT_NAME} plot.py)
python_install_on_site(${PROJECT_NAME} optimization.py)
