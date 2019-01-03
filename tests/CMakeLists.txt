cmake_minimum_required(VERSION 2.6)

include_directories("${PROJECT_SOURCE_DIR}/include")

FILE(GLOB_RECURSE HeaderFiles "${PROJECT_SOURCE_DIR}/include/spline/*.h")

add_executable(
				spline_tests ${HeaderFiles} Main.cpp
)

