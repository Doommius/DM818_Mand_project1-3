cmake_minimum_required(VERSION 3.5)
project(DM818)

#export CC = /opt/intel/bin/icc
set (CMAKE_CXX_STANDARD 11)
find_package(CUDA)
find_package(MPI REQUIRED)

include_directories(${MPI_INCLUDE_PATH})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

set(SOURCE_FILES_dgemm mandatory1/src/benchmark.cpp mandatory1/src/dgemm.cpp mandatory1/cuda/dev_array.h mandatory1/cuda/matrixmul.cu mandatory1/cuda/kernel.h mandatory1/cuda/Kernel.cu mandatory1/cuda/dgemm-cuda.cpp)
add_executable(dgemm ${SOURCE_FILES_dgemm})


set(SOURCE_FILES_testdgemm mandatory1/src/test.cpp mandatory1/src/dgemm.cpp)
add_executable(test-dgemm ${SOURCE_FILES_testdgemm})


set(SOURCE_FILES_mand2_serial mandatory2/src/serial/serial.cpp mandatory2/src/serial/common.cpp mandatory2/src/serial/common.h mandatory2/src/serial/grid.h mandatory2/src/serial/grid.cpp mandatory2/src/serial/serial.h)

add_executable(Particle_serial ${SOURCE_FILES_mand2_serial})

set(SOURCE_FILES_mand2_serial_opti mandatory2/src/serial/serial_opti.cpp mandatory2/src/serial/common.cpp  mandatory2/src/serial/grid.cpp mandatory2/src/serial/common.h  mandatory2/src/serial/serial.h mandatory2/src/serial/grid.h )
add_executable(Particle_serial_opti ${SOURCE_FILES_mand2_serial_opti})



set(SOURCE_FILES_mand2_serial_opti2 mandatory2/src/serial_opti/serial.cpp mandatory2/src/serial_opti/common.cpp mandatory2/src/serial_opti/common.h )
add_executable(Particle_serial_opti2 ${SOURCE_FILES_mand2_serial_opti2})