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
void unpackCBlock(double *C, unsigned int M, unsigned int N) {
    unsigned int c = 0;
    for (unsigned int n = 0; n < lda; n++) {
        for (unsigned int m = 0; m < lda; m++) {
            C[c++] = Cblock[m + n * lda];
        }
    }
}


//void printArray(double array[], int arraySize) {
//    for (int i = 0; i < arraySize; i++) {
//        std::cout << array[i] << ' ';
//    }
//    std::cout << '\n';
//}


//This is pretty what much Jacob did, except i wrote it with AVX2
//He said we should let the compiler take care of how to handle this. eg.
// __m256d c01x0-3 could be written as __m256d c23 = _mm_set_pd(0.0,0.0,0.0,0.0); and same goes for most of this function.
void core_4_4(double *A,double *B,double *C, unsigned int K){
    __m256d c01x0 = _mm256_setzero_pd();
    __m256d c23x0 = _mm256_setzero_pd();
    __m256d c01x1 = _mm256_setzero_pd();
    __m256d c23x1 = _mm256_setzero_pd();
    __m256d c01x2 = _mm256_setzero_pd();
    __m256d c23x2 = _mm256_setzero_pd();
    __m256d c01x3 = _mm256_setzero_pd();
    __m256d c23x3 = _mm256_setzero_pd();;


    for (unsigned int k = 0; k<K; k++){
        __m256d a01k = _mm256_load_pd(A+k*MR);
        __m256d a23k = _mm256_load_pd(A+2+k*MR);

        __m256d t1, t2;

        t1 = _mm256_set1_pd(B[k]);
        t2 = t1;
        t1 = _mm256_mul_pd(t1, a01k);
        c01x0 = _mm256_add_pd(c01x0, t1);
        t2 = _mm256_mul_pd(t2, a23k);
        c23x0 = _mm256_add_pd(c23x0, t2);

        t1 = _mm256_set1_pd(B[k+K]);
        t2 = t1;
        t1 = _mm256_mul_pd(t1, a01k);
        c01x1 = _mm256_add_pd(c01x1, t1);
        t2 = _mm256_mul_pd(t2, a23k);
        c23x1 = _mm256_add_pd(c23x1, t2);

        t1 = _mm256_set1_pd(B[k+2*K]);
        t2 = t1;
        t1 = _mm256_mul_pd(t1, a01k);
        c01x2 = _mm256_add_pd(c01x2, t1);
        t2 = _mm256_mul_pd(t2, a23k);
        c23x2 = _mm256_add_pd(c23x2, t2);

        t1 = _mm256_set1_pd(B[k+3*K]);
        t2 = t1;
        t1 = _mm256_mul_pd(t1, a01k);
        c01x3 = _mm256_add_pd(c01x3, t1);
        t2 = _mm256_mul_pd(t2, a23k);
        c23x3 = _mm256_add_pd(c23x3, t2);
    }
    _mm256_store_pd(C, c01x0);
    _mm256_store_pd(C+2, c23x0);

    _mm256_store_pd(C+MC, c01x0);
    _mm256_store_pd(C+2+MC, c23x0);

    _mm256_store_pd(C+2*MC, c01x0);
    _mm256_store_pd(C+2+2*MC, c23x0);

    _mm256_store_pd(C+3*MC, c01x0);
    _mm256_store_pd(C+2+3*MC, c23x0);

}

//maybe working, currently unsure.
void core_dyn(double *A,double *B,double *C, unsigned int M,unsigned int N,unsigned int K){
    for(unsigned int j = 0; j < N; j++){
        for (unsigned int i; i < M; i++){
            double cij = C[j*lda + i];
            for (int k = 0; k < K; ++k) {
                cij += A[k*lda + i] * B[j*lda + k];
            }
            C[j*lda + i] = cij;
        }
    }
}


void do_sub_block(double *A, double *B, double *C, unsigned int K, unsigned X, unsigned Y) {

};


void Prepare_block(double *C, unsigned int M, unsigned int N, unsigned int K) {
    double *B = Bblock;
    for (unsigned int n = 0; n < N; n += NR) {
        for (unsigned int m = 0; m < M; n += MR) {
            unsigned int Max_M = std::min(NR, M - m);
            unsigned int Max_N = std::min(MR, N - n);
            if (Max_M == MR && Max_N == NR){
                core_4_4(Ablock+m*K, B, Cblock+m,K);
            }else{
                core_dyn(Ablock+m*K,B,Cblock+m, Max_M, Max_N,K);
            }

        }
        B+= NR*K;
        unpackC(C+n*lda, M, std::min(NR, n-n));
    }
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
