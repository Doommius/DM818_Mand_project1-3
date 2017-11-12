//
// Created by jervelund on 10/26/17.
//

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <assert.h>
#include <time.h>
#define size 100   // Matrix size
#define cols size   // Matrix width
#define rows size   // Matrix height

void checkCUDAError(const char *msg)
{
    cudaError_t err = cudaGetLastError();
    if( cudaSuccess != err)
    {
        fprintf(stderr, "Cuda error: %s: %s.\n", msg, cudaGetErrorString( err) );
        exit(EXIT_FAILURE);
    }
}
__global__ void matrixMul( int *A, int *B, int *C)
{
    int bx = blockIdx.x; // Block index
    int tx = threadIdx.x; // Thread index
    int ts = blockDim.x; // number of threads
    // Declaration of the shared memory C element
    extern __shared__ int c_element_sum[];
    c_element_sum[tx] = A[tx+((bx/ts)*ts)] * B[(bx%ts)+(tx*ts)];

    //Block until all threads in the block have written their data to shared mem
    __syncthreads();

    int sum;
    for(int i=0; i<ts; i++){
        if(i==0){
            sum=c_element_sum[i];
        }
        else{
            sum+=c_element_sum[i];
        }
    }
    C[bx] = sum;

}


/////////////////////////////////////////////////////////
// Program main
/////////////////////////////////////////////////////////

int main(int argc, char** argv)
{
    //create timer.
    clock_t t1, t2;

    //start timer
    t1=clock();

    //allocate host memory for matrices
    unsigned int size_A = cols * rows;
    unsigned int mem_size_A = sizeof(int) * size_A;
    int* mA = (int*) malloc(mem_size_A);

    unsigned int size_B = cols * rows;
    unsigned int mem_size_B = sizeof(int) * size_B;
    int* mB = (int*) malloc(mem_size_B);

    unsigned int size_C = cols * rows;
    unsigned int mem_size_C = sizeof(int) * size_C;
    int* mC = (int*) malloc(mem_size_C);

    //initialize host memory
    for (int i = 0; i < size_A; ++i){
        mA[i] = 1;
        mB[i] = 1;
        mC[i] = 0;
    }

    // allocate device memory
    int* d_mA;
    int* d_mB;
    int* d_mC;
    cudaMalloc((void**) &d_mA, mem_size_A);
    cudaMalloc((void**) &d_mB, mem_size_B);
    cudaMalloc((void**) &d_mC, mem_size_C);

    //copy host memory to device (A and B)
    cudaMemcpy(d_mA, mA, mem_size_A, cudaMemcpyHostToDevice);
    cudaMemcpy(d_mB, mB, mem_size_B, cudaMemcpyHostToDevice);
    cudaMemcpy(d_mC, mC, mem_size_C, cudaMemcpyHostToDevice);

    // setup execution parameters
    int numThreadsPerBlock = cols;
    int numBlocks = (cols * rows);
    int sharedMemSize = numThreadsPerBlock * sizeof(int);

    dim3 dimGrid(numBlocks);
    dim3 dimBlock(numThreadsPerBlock);

    // execute the kernel
    matrixMul <<< dimGrid, dimBlock, sharedMemSize >>>(d_mA, d_mB, d_mC);

    //Block until device has completed
    cudaThreadSynchronize();

    // check if kernel execution generated an error
    // Check for any CUDA errors
    checkCUDAError("kernel invocation");

    //copy result from device to host
    cudaMemcpy(mC, d_mC, mem_size_C, cudaMemcpyDeviceToHost);

    // Check for any CUDA errors
    checkCUDAError("memcpy");

    //stop timer
    t2 = clock();

    //check results
    for (int i = 0; i < size_C; ++i){
        assert(mC[i] == cols);
    }

    //clean up memory
    free(mA);
    free(mB);
    free(mC);
    cudaFree(d_mA);
    cudaFree(d_mB);
    cudaFree(d_mC);

    printf("WITH CUDA - clocks: %d \n\n", t2-t1);

    //////////////////////////////
    ///////// CPU ONLY //////////
    /////////////////////////////

    //create timer.
    clock_t cpu_t1, cpu_t2;

    //start timer
    cpu_t1=clock();

    //allocate host memory for matrices
    unsigned int cpu_size_A = cols * rows;
    unsigned int cpu_mem_size_A = sizeof(int) * cpu_size_A;
    int* cpu_mA = (int*) malloc(cpu_mem_size_A);

    unsigned int cpu_size_B = cols * rows;
    unsigned int cpu_mem_size_B = sizeof(int) * cpu_size_B;
    int* cpu_mB = (int*) malloc(cpu_mem_size_B);

    unsigned int cpu_size_C = cols * rows;
    unsigned int cpu_mem_size_C = sizeof(int) * cpu_size_C;
    int* cpu_mC = (int*) malloc(cpu_mem_size_C);

    //initialize host memory
    for (int i = 0; i < cpu_size_A; ++i){
        cpu_mA[i] = 1;
        cpu_mB[i] = 1;
        cpu_mC[i] = 0;
    }

    int ts = cols;
    for(int bx=0; bx<(cols*rows);bx++){
        int sum = 0;
        for(int tx=0; tx<cols; tx++){
            sum += cpu_mA[tx+((bx/ts)*ts)] * cpu_mB[(bx%ts)+(tx*ts)];
        }
        cpu_mC[bx]=sum;
    }

    //stop timer
    cpu_t2 = clock();

    //check results
    for (int i = 0; i < cpu_size_C; ++i){
        assert(cpu_mC[i] == cols);
    }

    //clean up memory
    free(cpu_mA);
    free(cpu_mB);
    free(cpu_mC);

    printf("CPU ONLY - clocks: %d \n\n", cpu_t2-cpu_t1);

    return 0;
}