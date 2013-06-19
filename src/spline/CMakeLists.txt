cmake_minimum_required(VERSION 2.6)

add_definitions(-DSPLINE_DLLEXPORT)

file(
        GLOB_RECURSE
        source_files
        ./*
        ./API/*
)
 
add_library(
        spline
        SHARED
        ${source_files}
)
