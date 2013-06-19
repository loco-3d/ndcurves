cmake_minimum_required(VERSION 2.6)

include_directories("${PROJECT_SOURCE_DIR}/src/spline/API")

add_executable(
        spline_tests Main.cpp
)

target_link_libraries (spline_tests spline)
