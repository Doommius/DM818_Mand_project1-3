
#ifndef __COMM_H
#define __COMM_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

/* Function prototypes */
void spread_matrix(int processorCount, int matrixDimension, double *matrixA, double *matrixB);
void gather_matrix(/* PARAMETERS */);
int comms_split(/* PARAMETERS */);

#endif /* __COMM_H */
