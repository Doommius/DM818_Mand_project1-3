cmake_minimum_required(VERSION 3.5)
project(DM818)

set (CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx2")

set(SOURCE_FILES_dgemm mandatory1/src/benchmark.cpp mandatory1/src/dgemm.cpp)
add_executable(dgemm ${SOURCE_FILES_dgemm})


set(SOURCE_FILES_testdgemm mandatory1/src/test.cpp mandatory1/src/dgemm.cpp)
add_executable(test-dgemm ${SOURCE_FILES_testdgemm})