cmake_minimum_required(VERSION 3.5)
project(DM818)

#export CC = /opt/intel/bin/icc
set (CMAKE_CXX_STANDARD 11)
#find_package(CUDA)

find_package(MPI REQUIRED)
find_package(BLAS REQUIRED)

include_directories(${MPI_INCLUDE_PATH})

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

set(SOURCE_FILES_dgemm mandatory1/src/benchmark.cpp mandatory1/src/dgemm.cpp mandatory1/cuda/dev_array.h mandatory1/cuda/matrixmul.cu mandatory1/cuda/kernel.h mandatory1/cuda/Kernel.cu mandatory1/cuda/dgemm-cuda.cpp)
add_executable(mand_1_dgemm ${SOURCE_FILES_dgemm})


set(SOURCE_FILES_testdgemm mandatory1/src/test.cpp mandatory1/src/dgemm.cpp)
add_executable(mand_1_test-dgemm ${SOURCE_FILES_testdgemm})


set(SOURCE_FILES_mand2_serial mandatory2/src/serial/serial.cpp mandatory2/src/serial/common.cpp  mandatory2/src/serial/grid.cpp mandatory2/src/serial/common.h  )
add_executable(mand_2_serial ${SOURCE_FILES_mand2_serial})


set(SOURCE_FILES_mand2_serial_opti mandatory2/src/serial_opti/serial.cpp mandatory2/src/serial_opti/common.cpp mandatory2/src/serial_opti/common.h  mandatory2/src/serial_opti/grid.cpp mandatory2/src/serial_opti/grid.h)
add_executable(mand_2_serial_opti ${SOURCE_FILES_mand2_serial_opti})


set(SOURCE_FILES_mand2_openmp mandatory2/src/openmp/openmp.cpp mandatory2/src/openmp/common.cpp mandatory2/src/openmp/common.h)
add_executable(mand2_openmp ${SOURCE_FILES_mand2_openmp})

set(SOURCE_FILES_mand2_openmp_opti mandatory2/src/openmp_opti/openmp_opti.cpp mandatory2/src/openmp_opti/common.cpp mandatory2/src/openmp_opti/common.h mandatory2/src/openmp_opti/grid.h mandatory2/src/openmp_opti/grid.cpp)
add_executable(mand2_openmp_opti ${SOURCE_FILES_mand2_openmp_opti})



set(SOURCE_FILES_mand2_mpi mandatory2/src/mpi/mpi.cpp mandatory2/src/mpi/common.cpp mandatory2/src/mpi/common.h)
add_executable(mand2_mpi ${SOURCE_FILES_mand2_mpi})

set(SOURCE_FILES_mand2_mpi_opti mandatory2/src/mpi_opti/mpi.cpp mandatory2/src/mpi_opti/common.cpp mandatory2/src/mpi_opti/common.h mandatory2/src/mpi_opti/grid.cpp mandatory2/src/mpi_opti/grid.h)
add_executable(mand2_mpi_opti ${SOURCE_FILES_mand2_mpi_opti})

target_link_libraries(mand2_mpi ${MPI_LIBRARIES})

target_link_libraries(mand2_mpi_opti ${MPI_LIBRARIES})


set(SOURCE_FILES_mand3 mandatory3/src/main.c mandatory3/src/matrix.c mandatory3/src/matrix.h)
add_executable(mand3 ${SOURCE_FILES_mand3})

target_link_libraries(mand3 ${MPI_LIBRARIES} ${BLAS_LIBRARIES} m )



#flag for openmp -fopenmp