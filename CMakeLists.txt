cmake_minimum_required(VERSION 3.5)
project(DM818)

set (CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES src/benchmark.cpp src/dgemm-blocked.cpp src/dgemm-naive.cpp src/dgemm-blas.cpp src/dgemm.cpp)
set(SOURCE_FILES_naive src/benchmark.cpp src/dgemm-naive.cpp)
set(SOURCE_FILES_blocked src/benchmark.cpp src/dgemm-blocked.cpp)
set(SOURCE_FILES_blas src/benchmark.cpp src/dgemm-blas.cpp)
add_executable(DM818 ${SOURCE_FILES})


add_executable(benchmarknaive ${SOURCE_FILES_naive})
add_executable(benchmarkblocked ${SOURCE_FILES_blocked})
add_executable(benchmarkblas ${SOURCE_FILES_blas})
