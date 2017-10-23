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

//debugging
void printMatrix2(double *matrix, int MatrixSize) {
    std::cout << ("Starting Matrix") << '\n';
    for (int i = 0; i < MatrixSize * MatrixSize; i++) {
        std::cout << matrix[i] << ' ';
        if ((i + 1) % MatrixSize == 0) {
            printf("\n");
        }
    }
    std::cout << ("Ending matrix") << '\n';
}

void printdouble(char *text, double *X, int arraySize) {
    std::cout << text << '\n';
    for (int i = 0; i < arraySize; i++) {
        double a = X[i];

        std::cout << a << ' ';
    }
    std::cout << '\n';
}


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


//TODO currently broken, currently fixed
void unpackCBlock(double *C, unsigned int M, unsigned int N) {
//    printMatrix2(Cblock,4 );
    for (unsigned int n = 0; n < NR; n++) {
        for (unsigned int i = 0; i < M; i++) {
//            int lookingat = i + n * MC;    ////Debugging value used when debugging with gdb
//            int storeingat = i + n * lda;  ////Debugging value used when debugging with gdb
//            double x = Cblock[lookingat];  ////Debugging value used when debugging with gdb
            C[i + n * lda] = Cblock[i + n * MC];
        }
    }
}


//This is pretty what much Jacob did, except i wrote it with AVX2
//He said we should let the compiler take care of how to handle this. eg.
// __m256d c01x0-3 could be written as __m256d c23 = _mm_set_pd(0.0,0.0,0.0,0.0); and same goes for most of this function.
void core_4_4(double *A, double *B, double *C, unsigned int K) {
//     std::cout << "K is = " << K << '\n';
    __m256d c0 = _mm256_setzero_pd();
    __m256d c1 = _mm256_setzero_pd();
    __m256d c2  = _mm256_setzero_pd();
    __m256d c3  = _mm256_setzero_pd();

    for (unsigned int k = 0; k < K - 1; k += 2) {
//        printdouble(A+k*MR,4);
//        printdouble(A+(2+k)*MR,4);

        __m256d a1 = _mm256_load_pd(A + k * MR);
//        printdouble("A0123k0",A + (k) * MR, 4);

        __m256d a2 = _mm256_load_pd(A + (k + 1) * MR);

//        printdouble("A0123k1",A+ (k+1) * MR, 4);

        __m256d b00  = _mm256_set1_pd(B[k]);
        __m256d b01  = _mm256_set1_pd(B[k + K]);
        __m256d b02  = _mm256_set1_pd(B[k + K * 2]);
        __m256d b03  = _mm256_set1_pd(B[k + K * 3]);

        __m256d b10  = _mm256_set1_pd(B[k + 1]);
        __m256d b11  = _mm256_set1_pd(B[k + K + 1]);
        __m256d b12  = _mm256_set1_pd(B[k + K * 2 + 1]);
        __m256d b13  = _mm256_set1_pd(B[k + K * 3 + 1]);

        //do left
        c0 = _mm256_add_pd(c0, _mm256_mul_pd(b00, a1));
        c1 = _mm256_add_pd(c1, _mm256_mul_pd(b01, a1));
        c2 = _mm256_add_pd(c2, _mm256_mul_pd(b02, a1));
        c3 = _mm256_add_pd(c3, _mm256_mul_pd(b03, a1));

        c0 = _mm256_add_pd(c0, _mm256_mul_pd(b10, a2));
        c1 = _mm256_add_pd(c1, _mm256_mul_pd(b11, a2));
        c2 = _mm256_add_pd(c2, _mm256_mul_pd(b12, a2));
        c3 = _mm256_add_pd(c3, _mm256_mul_pd(b13, a2));

    }
    //edge case.
    if (K % 2 == 1) {
        unsigned int k = K - 1;
        __m256d a0 = _mm256_load_pd(A + k * MR);

        __m256d b00 = _mm256_set1_pd(B[k]);
        __m256d b01 = _mm256_set1_pd(B[k + K]);
        __m256d b02 = _mm256_set1_pd(B[k + K * 2]);
        __m256d b03 = _mm256_set1_pd(B[k + K * 3]);

        c0 = _mm256_add_pd(c0, _mm256_mul_pd(b00, a0));
        c1 = _mm256_add_pd(c1, _mm256_mul_pd(b01, a0));
        c2 = _mm256_add_pd(c2, _mm256_mul_pd(b02, a0));
        c3 = _mm256_add_pd(c3, _mm256_mul_pd(b03, a0));
    }

    //TODO something is broken here. its been fixed
    ////Debugging value used when debugging with gdb
//    double *tmp = (double *) _mm_malloc(NR * sizeof(double), 32);
//
//    _mm256_store_pd(tmp, c0123x0);
//
//    _mm256_store_pd(tmp, c0123x1);
//
//    _mm256_store_pd(tmp, c0123x2);
//
//    _mm256_store_pd(tmp, c0123x3);

//    _mm_free(tmp);

    _mm256_store_pd(C, c0);
    _mm256_store_pd(C + MC, c1);
    _mm256_store_pd(C + 2 * MC, c2);
    _mm256_store_pd(C + 3 * MC, c3);


}

//works.
void core_dyn(double *A, double *B, double *C, unsigned int M, unsigned int N, unsigned int K) {
    for (unsigned int j = 0; j < N; j++) {
        for (unsigned int i = 0; i < M; i++) {
            double cij = C[j * lda + i];
            for (int k = 0; k < K; ++k) {
                cij += A[k * lda + i] * B[j * lda + k];
            }
            C[j * MC + i] += cij;
        }
    }
}

void Do_block(double *C, unsigned int M, unsigned int N, unsigned int K) {
    double *B = Bblock;
    for (unsigned int n = 0; n < N; n += NR) {
        for (unsigned int m = 0; m < M; m += MR) {
            unsigned int Max_M = std::min(NR, M - m);
            unsigned int Max_N = std::min(MR, N - n);
            if (Max_M == MR && Max_N == NR) {
                core_4_4(Ablock + m * K, B, Cblock + m, K);
            } else {
                core_dyn(Ablock + m, B, Cblock +m, M, N, K);
            }
        }
        B += NR * K;
        unpackCBlock(C + n * lda, M, std::min(NR, n - n));
    }
}

void square_dgemm(int M, double *A, double *B, double *C) {
    lda = M;
    Ablock = (double *) _mm_malloc(MC * KC * sizeof(double), 32); //128*128
    Bblock = (double *) malloc(lda * KC * sizeof(double)); // M * 128
    Cblock = (double *) _mm_malloc(MC * NR * sizeof(double), 32); //128*4

    for (unsigned int k = 0; k < lda; k += KC) {
        packBBlock(B + k, std::min(KC, lda - k));
        for (unsigned int i = 0; i < lda; i += MC) {
            packAblock(A + i + k * lda, std::min(MC, lda - i), std::min(KC, lda - k));
            Do_block(C + i, std::min(MC, lda - i), lda, std::min(KC, lda - k));
        }
    }
    _mm_free(Ablock);
    free(Bblock);
    _mm_free(Cblock);
}
