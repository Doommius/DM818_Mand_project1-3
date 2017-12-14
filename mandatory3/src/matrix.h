
#ifndef MATRIX_H
#define MATRIX_H
#include <stdlib.h>
#include <mpi.h>

/* Function prototypes */
void make_full_matrices(int matrixDimensions, double *matrixA, double *matrixB);
void allocate_matrices(int matrixDimensions, double *matrixA, double *matrixB,double *matrixC);
void reset_matrices(int matrixDimensions, double *matrixA, double *matrixB,double *matrixC);
void destroy_matrices(double *matrixA, double *matrixB,double *matrixC  );
void matrix_mult(int N, double* A, double* B, double* C);
void setzero(int matrixDimensions, double *matrixC);
void sumMatrices(void *in, void *inout, int *length, MPI_Datatype *type);


#endif /* MATRIX_H */
