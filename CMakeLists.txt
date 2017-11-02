cmake_minimum_required(VERSION 3.5)
project(DM818)

#export CC = /opt/intel/bin/icc
set (CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx2 -O3")

set(SOURCE_FILES_dgemm mandatory1/src/benchmark.cpp mandatory1/src/dgemm.cpp)
add_executable(dgemm ${SOURCE_FILES_dgemm})


set(SOURCE_FILES_testdgemm mandatory1/src/test.cpp mandatory1/src/dgemm.cpp)
add_executable(test-dgemm ${SOURCE_FILES_testdgemm})




set(SOURCE_FILES_mand2_serial mandatory2/src/serial.cpp /mandatory2/src/common.cpp /mandatory2/src/common.h)
add_executable(Particle_serial ${SOURCE_FILES_mand2_serial})


set(SOURCE_FILES_mand2_serial_opti mandatory2/src/serial_opti.cpp /mandatory2/src/common.cpp /mandatory2/src/common.h mandatory2/src/grid.cpp mandatory2/src/grid.h)
add_executable(Particle_serial_opti ${SOURCE_FILES_mand2_serial_opti})




