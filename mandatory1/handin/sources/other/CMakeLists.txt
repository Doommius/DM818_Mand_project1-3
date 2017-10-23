cmake_minimum_required(VERSION 3.5)
project(DM818)

set (CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx2")


set(SOURCE_FILES benchmark.cpp dgemm-blocked.cpp dgemm-naive.cpp dgemm-blas.cpp dgemm.cpp test.cpp)
add_executable(DM818 ${SOURCE_FILES})


set(SOURCE_FILES_naive benchmark.cpp dgemm-naive.cpp)
add_executable(benchmarknaive ${SOURCE_FILES_naive})


set(SOURCE_FILES_blocked benchmark.cpp dgemm-blocked.cpp)
add_executable(benchmarkblocked ${SOURCE_FILES_blocked})


set(SOURCE_FILES_blas benchmark.cpp dgemm-blas.cpp)
add_executable(benchmarkblas ${SOURCE_FILES_blas})



set(SOURCE_FILES_testmnaive test.cpp dgemm-naive.cpp)
add_executable(Naivetest ${SOURCE_FILES_testmnaive})


set(SOURCE_FILES_dgemm benchmark.cpp dgemm.cpp)
add_executable(dgemm ${SOURCE_FILES_dgemm})


set(SOURCE_FILES_testdgemm test.cpp dgemm.cpp)
add_executable(test-dgemm ${SOURCE_FILES_testdgemm})