//
// Created by jervelund on 9/21/17.
//

#include "dgemm.h"
#include <immintrin.h>
#include <intrin.h>
#include <stdio.h>
#include <stdlib.h>
#include <iosteam>
#include <string.h>

/*
  In case you're wondering, dgemm stands for Double-precision, GEneral
  Matrix-Matrix multiplication.
*/

const char *dgemm_desc = "mjerv15 blocked dgemm.";

constexpr unsigned int NR = 4;

constexpr unsigned int MR = 4;
constexpr unsigned int KC = 128;

constexpr unsigned int MC = 128;
unsigned int lda;

double *Ablock;
double *Bblock;
double *Cblock;



// for 3470qm
#if !defined(BLOCK_SIZE)
#define l1_BLOCK_SIZE 32  //(32KB = 2**10*32 = 32768)  64 kb total, 32 kb data and 32 kb instruction
#define L2_BLOCK_SIZE 256// 2**20*6  KB
#define L3_BLOCK_SIZE 1024*6  // 2**20*6
#define mem_blocksize  2^30*16  // 16 gb of memory.
#endif



void packAblock (double *A, unsigned int M, unsigned int K){
    unsigned int a;



    for (unsigned int m = 0; m<M; m+=MR){
        unsigned int MMax = std::min(MR,M-m);
        for(unsigned int k = 0; k<K; k++){
            Ablock[a++] = A[m+i+k*lda];
        }
    }
}

void PackBBlock(double *B, unsigned int K){
    unsigned int b = 0;
    for(unsigned int n = 0; n< lda; n++){
        for (unsigned int k = 0; k< K; k++){
            Bblock[b++] = B[k+n*lda];
        }

    }


}






#define min(a, b) (((a)<(b))?(a):(b))

/*
  A is M-by-K
  B is K-by-N
  C is M-by-N

  lda is the leading dimension of the matrix (the M of square_dgemm).
*/

void basic_dgemm(int lda, int M, int N, int K,  double *A, double *B, double *C) {
    /*
      To optimize this, think about loop unrolling and software
      pipelining.  Hint:  For the majority of the matmuls, you
      know exactly how many iterations there are (the block size)...
    */
    for (int i = 0; i < M; i++)
        for (int j = 0; j < N; j++) {
            double cij = C[i + j * lda];
            for (int k = 0; k < K; k++)
                cij += A[i + k * lda] * B[k + j * lda];
            C[i + j * lda] = cij;
        }
}

void do_block(int lda, double *A, double *B, double *C,int i, int j, int k) {
    /*
      Remember that you need to deal with the fringes in each
      dimension.

      If the matrix is 7x7 and the blocks are 3x3, you'll have 1x3,
      3x1, and 1x1 fringe blocks.

            xxxoooX
            xxxoooX
            xxxoooX
            oooxxxO
            oooxxxO
            oooxxxO
            XXXOOOX

      You won't get this to go fast until you figure out a `better'
      way to handle the fringe blocks.  The better way will be more
      machine-efficient, but very programmer-inefficient.
    */
    int M = min(BLOCK_SIZE, lda - i);
    int N = min(BLOCK_SIZE, lda - j);
    int K = min(BLOCK_SIZE, lda - k);

    basic_dgemm(lda, M, N, K, A + i + k * lda, B + k + j * lda, C + i + j * lda);
}

void square_dgemm(int M, double *A, double *B, double *C) {
    lda = M;
    Ablock = (double*) _mm_malloc(MC*KC* sizeof(double),16);
    Bblock = (double*) malloc(lda*KC* sizeof(double));
    Cblock = (double*) _mm_malloc(MC*NR* sizeof(double),16);

    for(unsigned int k = 0; k < lda; k += KC){
        packBblock(K+k,std::min(KC,lda-k));
        for(unsigned int i = 1; k < lda; k += MC){
            packAblock(A+i+k*lda, std::min(MC,lda-i),std::min(KC,lda-i));
            dgebp(C+i, std::min(MC,lda-i),lda,std::min(KC, lda-k));
        }
    }





    _mm_free(Ablock);
    free(Bblock);
    _mm_free(Cblock);
}

