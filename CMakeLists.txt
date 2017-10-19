cmake_minimum_required(VERSION 3.5)
project(DM818)

set (CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx2")


set(SOURCE_FILES src/benchmark.cpp src/dgemm-blocked.cpp src/dgemm-naive.cpp src/dgemm-blas.cpp src/dgemm-mjerv.cpp src/dgemm-mjerv2.cpp src/test.cpp)
set(SOURCE_FILES_naive src/benchmark.cpp src/dgemm-naive.cpp)
set(SOURCE_FILES_blocked src/benchmark.cpp src/dgemm-blocked.cpp)
set(SOURCE_FILES_blas src/benchmark.cpp src/dgemm-blas.cpp src/dgemm-mjervbackup.cpp)




add_executable(DM818 ${SOURCE_FILES})
add_executable(benchmarknaive ${SOURCE_FILES_naive})
add_executable(benchmarkblocked ${SOURCE_FILES_blocked})
add_executable(benchmarkblas ${SOURCE_FILES_blas})



set(SOURCE_FILES_testmnaive src/test.cpp src/dgemm-naive.cpp)
add_executable(Naivetest ${SOURCE_FILES_testmnaive})


set(SOURCE_FILES_mjerv15 src/benchmark.cpp src/dgemm-mjerv.cpp)
add_executable(mjerv ${SOURCE_FILES_mjerv15})


set(SOURCE_FILES_testmjerv src/test.cpp src/dgemm-mjerv.cpp)
add_executable(mjervtest ${SOURCE_FILES_testmjerv})