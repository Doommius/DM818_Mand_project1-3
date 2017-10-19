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
        p[i] = i+1;
}
void filltest0(double *p, int n) {
    for (int i = 0; i < n; i++)
        p[i] = 0;
}

void printMatrix(double *matrix, int MatrixSize) {
    std::cout << ("Starting Matrix") << '\n';
    for (int i = 0; i < MatrixSize*MatrixSize; i++) {
        std::cout << matrix[i] << ' ';
        if ((i+1) % MatrixSize == 0){
            printf("\n");
        }
    }
    std::cout << ("Ending matrix") << '\n';



}


void naive_square_dgemm( int n, double *A, double *B, double *C )
{
    for( int i = 0; i < n; i++ )
        for( int j = 0; j < n; j++ )
        {
            double cij = C[i+j*n];
            for( int k = 0; k < n; k++ )
                cij += A[i+k*n] * B[k+j*n];
            C[i+j*n] = cij;
        }
}



int main(int argc, char **argv) {

    srand(time(NULL));
    printf("Test Program running \n");
//    printf("Description:\t%s\n\n", dgemm_desc);

    //
    // These sizes should highlight performance dips at multiples of certain
    // powers-of-two
    //
    int testworks_sizes[] = {4};
    int test_sizes[] = {3};
//
//    int test_sizes[] = {2};

    for (int isize = 0; isize < sizeof(test_sizes) / sizeof(test_sizes[0]); isize++) {
        int n = test_sizes[isize];
        int n1 = n;

        double *A = (double *) malloc(n * n * sizeof(double));
        double *B = (double *) malloc(n * n * sizeof(double));
        double *C = (double *) malloc(n * n * sizeof(double));
        double *A1 = (double *) malloc(n * n * sizeof(double));
        double *B1 = (double *) malloc(n * n * sizeof(double));
        double *C1 = (double *) malloc(n * n * sizeof(double));

        filltest(A, n * n);
        filltest(B, n * n);
        filltest0(C, n * n);
        filltest(A1, n * n);
        filltest(B1, n * n);
        filltest0(C1, n * n);
        printf("Matrix A \n");
        printMatrix(A, n);
        printf("Matrix B \n");
        printMatrix(B, n);
        square_dgemm(n, A, B, C);
        printf("Matrix C (result) \n");
        printMatrix(C, n);
        naive_square_dgemm(n1,A1,B1,C1);
        printMatrix(C1,n1);
        free(A);
        free(B);
        free(C);

    }
    return 0;
}