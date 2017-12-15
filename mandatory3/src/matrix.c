
#include <stdlib.h>
#include "matrix.h"
#include <cblas.h> //assumes general CBLAS interface


//const char* dgemm_desc = "BLAS dgemm.";

#define DGEMM dgemm_

extern void DGEMM(char *, char *, int *, int *, int *, double *, double *, int *, double *, int *, double *, double *, int *);

/* Make P(0,0,0)'s full size matrices */
void make_full_matrices(int matrixDimensions, double *matrixA, double *matrixB) {
    for (int i = 0; i < matrixDimensions * matrixDimensions; i++) {
        matrixA[i] = 2 * drand48() - 1;
        matrixB[i] = 2 * drand48() - 1;
    }
}

/* Allocate matrices */
void allocate_matrices(int matrixDimensions, double *matrixA, double *matrixB, double *matrixC) {
    matrixA = malloc(matrixDimensions * matrixDimensions * sizeof(double));
    matrixB = malloc(matrixDimensions * matrixDimensions * sizeof(double));
    matrixC = malloc(matrixDimensions * matrixDimensions * sizeof(double));
}

void setzero(int matrixDimensions, double *matrixC) {
    for (int i = 0; i < matrixDimensions * matrixDimensions; i++) {
        matrixC[i] = 0.0;
    }
}

/* Reset the matrices for new use */
void reset_matrices(int matrixDimensions, double *matrixA, double *matrixB, double *matrixC) {
    for (int i = 0; i < matrixDimensions * matrixDimensions; i++) {
        matrixA[i] = 2 * drand48() - 1;
        matrixB[i] = 2 * drand48() - 1;
        matrixC[i] = 0.0;
    }
}

/* Free the memory for matrices */
void destroy_matrices(double *matrixA, double *matrixB, double *matrixC) {
    free(matrixA);
    free(matrixB);
    free(matrixC);
}

/* Matrix multiplication */
void matrix_mult(int N, double *A, double *B, double *C) {



    {
        char TRANSA = 'N';
        char TRANSB = 'N';
        int M = N;
        int K = N;
        double ALPHA = 1.;
        double BETA = 1.;
        int LDA = N;
        int LDB = N;
        int LDC = N;
        DGEMM(&TRANSA, &TRANSB, &M, &N, &K, &ALPHA, A, &LDA, B, &LDB, &BETA, C, &LDC);
    }


}


void naive_dgemm( int n, double *A, double *B, double *C )
{
//    for (int i = 0; i < N; i++) {
//        for (int j = 0; j < N; j++) {
//            double cij = C[i + j * N];
//            for (int k = 0; k < N; k++) {
//                cij += A[i + k * N] * B[k + j * N];
//            }
//            C[i + j * N] = cij;
//        }
//    }
//
//    return;
}

void sumMatrices(void *in, void *inout, int *length, MPI_Datatype *type) {
    double *a = (double *) in;
    double *b = (double *) inout;
    for (int i = 0; i < *length; i++) {
        b[i] += a[i];
    }
}

