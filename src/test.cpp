//
// Created by Mark on 14-Oct-17.
//
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <malloc.h>
#include <iostream>
#include <time.h>       /* time */


//
//  Your function must have the following signature:
//
extern const char* dgemm_desc;
extern void square_dgemm( int M, double *A, double *B, double *C );

//


void fill(double *p, int n) {
    for (int i = 0; i < n; i++)
        p[i] = 2 * drand48() - 1;
}

void filltest(double *p, int n) {
    for (int i = 0; i < n; i++)
        p[i] = 5;
}
void filltest0(double *p, int n) {
    for (int i = 0; i < n; i++)
        p[i] = 0;
}

void printMatrix(double *matrix, int MatrixSize) {
    for (int i = 0; i < MatrixSize; i++) {
        for (int i = 0; i < MatrixSize; i++) {
            std::cout << matrix[i] << ' ';
        }
        printf("\n");
    }
    std::cout << '\n';
}


int main(int argc, char **argv) {

    srand (time(NULL));
printf("Test Program running \n");
//    printf("Description:\t%s\n\n", dgemm_desc);

    //
    // These sizes should highlight performance dips at multiples of certain
    // powers-of-two
    //
    int test_sizes[] = {2};
//
//    int test_sizes[] = {2};

    for (int isize = 0; isize < sizeof(test_sizes) / sizeof(test_sizes[0]); isize++) {
        int n = test_sizes[isize];

        double *A = (double *) malloc(n * n * sizeof(double));
        double *B = (double *) malloc(n * n * sizeof(double));
        double *C = (double *) malloc(n * n * sizeof(double));

        filltest(A, n * n);
        filltest(B, n * n);
        filltest0(C, n * n);
        printf("Matrix A \n");
        printMatrix(A, n);
        printf("Matrix B \n");
        printMatrix(B, n);
        square_dgemm(n, A, B, C);
        printf("Matrix C (result) \n");
        printMatrix(C, n);
        return 0;
    }
}
