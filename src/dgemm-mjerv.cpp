//
// Created by jervelund on 9/21/17.
//


#include <immintrin.h>
#include <stdio.h>
#include <algorithm>
#include <malloc.h>
#include <iostream>

/*
  In case you're wondering, dgemm stands for Double-precision, GEneral  Matrix-Matrix multiplication.
*/

const char *dgemm_desc = "mjerv15 blocked dgemm.";

//#define min(a,b) (((a)<(b))?(a):(b))

unsigned int NR = 4;
unsigned int MR = 4;

unsigned int KC = 128;
unsigned int MC = 128;
unsigned int lda; //size lda*lda of matrix

double *Ablock;
double *Bblock;
double *Cblock;

//              From Lecture by Jacob
//		 				                      N
//                                +-------------------------+
//                                |                         |
//                                |                         |
//                                |                         |
//                                |                         |
//                                |                         |
//                              K |            B            |
//                                |                         |
//                                |                         |
//                                |                         |
//                                |                         |
//                                |                         |
//                                |                         |
// 				                  +-------------------------+
//		        K
//  +-------------------------+   +-------------------------+
//  |                         |   |                         |
//  |                         |   |                         |
//  |                         |   |                         |
//  |                         |   |                         |
//  |                         |   |                         |
//  |           A             |   |                         |
// M|                         |   |            C            |
//  |                         |   |                         |
//  |                         |   |                         |
//  |                         |   |                         |
//  |                         |   |                         |
//  |                         |   |                         |
//  +-------------------------+   +-------------------------+


void packAblock(double *A, unsigned int M, unsigned int K) {
    unsigned int a = 0;
    for (unsigned int m = 0; m < M; m += MR) {
        unsigned int MMax = std::min(MR, M - m);
        for (unsigned int k = 0; k < K; k++) {
            for (unsigned int i = 0; i < MMax; i++) {
                Ablock[a++] = A[m + i + k * lda];
            }

        }
    }
}

void packBBlock(double *B, unsigned int K) {
    unsigned int b = 0;
    for (unsigned int n = 0; n < lda; n++) {
        for (unsigned int k = 0; k < K; k++) {
            Bblock[b++] = B[k + n * lda];
        }

    }


}

void printArray(double array[], int arraySize) {
    for (int i = 0; i < arraySize; i++) {
        std::cout << array[i] << ' ';
    }
    std::cout << '\n';
}

void do_sub_block(double *A, double *B, double *C, unsigned int K, unsigned X, unsigned Y) {

};


void Prepare_block(double *C, unsigned int M, unsigned int N, unsigned int K) {
    double *B = Bblock;
    for (unsigned int n = 0; n < N; n += NR) {
        for (unsigned int m = 0; m < M; n += MR) {
            unsigned int Max_M = std::min(NR, M - m);
            unsigned int Max_N = std::min(MR, N - n);

            do_sub_block(Ablock, Bblock, C, K, Max_M, Max_N);
        }
        //unpackC(C+n*lda, M, std::min(NR, n-n));
    }
}

void core_4_4(double *A,double *B,double *C, unsigned int k){
    __m256d Vec10 = _mm256_set_pd(0.0,0.0,0.0,0.0);
    __m256d Vec11 = _mm256_set_pd(0.0,0.0,0.0,0.0);




    }
void square_dgemm(int M, double *A, double *B, double *C) {
    lda = M;
    Ablock = (double *) _mm_malloc(MC * KC * sizeof(double), 16); //128*128
    Bblock = (double *) malloc(lda * KC * sizeof(double)); // M * 128
    Cblock = (double *) _mm_malloc(MC * NR * sizeof(double), 16); //128*4

    for (unsigned int k = 0; k < lda; k += KC) {
        packBBlock(B + k, std::min(KC, lda - k));
        for (unsigned int i = 0; i < lda; i += MC) {
            packAblock(A + i + k * lda, std::min(MC, lda - i), std::min(KC, lda - k));
            Prepare_block(C + i, std::min(MC, lda - i), lda, std::min(KC, lda - k));
        }
    }
    _mm_free(Ablock);
    free(Bblock);
    _mm_free(Cblock);
}
