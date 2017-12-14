#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include "matrix.h"


int rank;
int myrank;
int maxrank;


MPI_Op matrixSum;


int coordinates[3];
MPI_Comm iComm, jComm, kComm, ijComm;

double *receivedMatrixA;
double *receivedMatrixB;
double *BlockMatrixC;

int blockLength;


/* Print a header for results output */
void results_header() {
    printf("Dims  No. Proc.  Avg. RT / Dev. (Eff.)\n");
}

/* Print the stats for 1 run */
void write_result(int full_dim, int procs, double rt, double dev, double eff) {
    printf("%-5i %-10i %-5.5f / %-5.5f (%-5.5f)\n", full_dim, procs, rt, dev, eff);
}

/* Average and standard deviation */
double average(int count, const double *list, double *dev) {
    int i;
    double sum = 0.0, avg;

    for (i = 0; i < count; i++) {
        sum += list[i];
    }

    avg = sum / (double) count;

    if (dev != 0) {
        sum = 0.0;
        for (i = 0; i < count; i++) {
            sum += (list[i] - avg) * (list[i] - avg);
        }

        *dev = sqrt(sum / (double) count);
    }

    return avg;
}

/* Divide and send/recieve the matrix before calculation */
void spread_matrix(int processorCount, int matrixDimension, double *matrixA, double *matrixB) {
    double *preparedMatrixA = NULL;
    double *preparedMatrixB = NULL;

    int sendCount[processorCount];
    int displacements[processorCount];

    int length = (int) cbrt(processorCount);
    int blockLength = matrixDimension / length;

    receivedMatrixA = (double *) malloc(sizeof(double) * blockLength * blockLength);
    receivedMatrixB = (double *) malloc(sizeof(double) * blockLength * blockLength);

    if (coordinates[2] == 0) { // Only perform initial distribution to k = 0
        if (rank == 0) {
            preparedMatrixA = (double *) malloc(sizeof(double) * matrixDimension * matrixDimension);
            preparedMatrixB = (double *) malloc(sizeof(double) * matrixDimension * matrixDimension);

            for (int i = 0; i < length; i++) {
                for (int j = 0; j < length; j++) {
                    for (int k = 0; k < blockLength; k++) {
                        // Offset into the prepared matrix
                        int offsetPrepared = i * length * (blockLength * blockLength) +
                                             j * (blockLength * blockLength) +
                                             k * blockLength;

                        // The start of the matrix
                        int offsetMatrix = j * matrixDimension * blockLength +
                                           i * blockLength +
                                           k * matrixDimension;

                        // Copy them into the prepared matrices
                        memcpy(&preparedMatrixA[offsetPrepared], &matrixA[offsetMatrix], sizeof(double) * blockLength);
                        memcpy(&preparedMatrixB[offsetPrepared], &matrixB[offsetMatrix], sizeof(double) * blockLength);
                    }
                }
            }
        }

        for (int i = 0; i < processorCount; i++) {
            sendCount[i] = blockLength * blockLength;
            displacements[i] = i * blockLength * blockLength;
        }

        // Perform the communication:
        //TODO Fix error.
        MPI_Scatterv(preparedMatrixA, sendCount, displacements, MPI_DOUBLE, receivedMatrixA,
                     blockLength * blockLength, MPI_DOUBLE, 0, ijComm);
        MPI_Scatterv(preparedMatrixB, sendCount, displacements, MPI_DOUBLE, receivedMatrixB,
                     blockLength * blockLength, MPI_DOUBLE, 0, ijComm);

        if (rank == 0) {
            free(preparedMatrixA);
            free(preparedMatrixB);
        }

    }
}

void distribute() {
    // Distribute matrix A
    if (coordinates[2] == 0 && coordinates[1] != 0) {
        MPI_Send(receivedMatrixA, blockLength * blockLength, MPI_DOUBLE, coordinates[1], 0, kComm);
    } else if (coordinates[1] == coordinates[2] && coordinates[1] != 0) {
        MPI_Recv(receivedMatrixA, blockLength * blockLength, MPI_DOUBLE, 0, 0, kComm, MPI_STATUS_IGNORE);
    } /* else do nothing */

    if (coordinates[2] == 0 && coordinates[0] != 0) {
        MPI_Send(receivedMatrixB, blockLength * blockLength, MPI_DOUBLE, coordinates[0], 0, kComm);
    } else if (coordinates[0] == coordinates[2] && coordinates[2] != 0) {
        MPI_Recv(receivedMatrixB, blockLength * blockLength, MPI_DOUBLE, 0, 0, kComm, MPI_STATUS_IGNORE);
    } /* else do nothing */
}

/**
 * Step C
 * Broadcast elements along j dimension for A and i dimension for B.
 */
void broadcast() {
    MPI_Bcast(receivedMatrixA, blockLength * blockLength, MPI_DOUBLE, coordinates[2], jComm);
    MPI_Bcast(receivedMatrixB, blockLength * blockLength, MPI_DOUBLE, coordinates[2], iComm);
}

void gather_matrix(double *matrixC) {
    MPI_Reduce(matrixC, BlockMatrixC, blockLength * blockLength, MPI_DOUBLE, matrixSum, 0, kComm);
}

/*
 * Step D
 * Multiply local process' A and B together.
 * Then reduce all k>0 to k=0
 */
void multiplyAndReduce() {
    double *matrixC = (double *) malloc(sizeof(double) * blockLength * blockLength);
    memset(matrixC, 0, sizeof(double) * blockLength * blockLength);
    matrix_mult(blockLength, receivedMatrixA, receivedMatrixB, matrixC);
    gather_matrix(matrixC);
    free(matrixC);
}

void initMPI(int argc, char **argv) {


    MPI_Comm gridCommunicator;

    MPI_Init(&argc, &argv); //TODO This causes the program to crash...

    MPI_Comm_size(MPI_COMM_WORLD, &maxrank);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int pEachDimension = (int) cbrt(maxrank);
    int dimensions[3] = {pEachDimension, pEachDimension, pEachDimension};
    int periods[3] = {0, 0, 0};

    if ((pEachDimension * pEachDimension * pEachDimension) != maxrank) {
        if (rank == 0) printf("Error! Number of processes is not a 3dim cube! Exiting..\n");
        MPI_Finalize();
        exit(-1);
    }

    MPI_Cart_create(MPI_COMM_WORLD, 3, dimensions, periods, 0, &gridCommunicator);
    MPI_Comm_rank(gridCommunicator, &rank);
    MPI_Cart_coords(gridCommunicator, rank, 3, coordinates);

    int iDimensions[3] = {1, 0, 0};
    int jDimensions[3] = {0, 1, 0};
    int kDimensions[3] = {0, 0, 1};
    int ijDimensions[3] = {1, 1, 0};

    MPI_Cart_sub(gridCommunicator, iDimensions, &iComm);
    MPI_Cart_sub(gridCommunicator, jDimensions, &jComm);
    MPI_Cart_sub(gridCommunicator, kDimensions, &kComm);
    MPI_Cart_sub(gridCommunicator, ijDimensions, &ijComm);

    MPI_Op_create(sumMatrices, 1, &matrixSum);
}


int main(int argc, char **argv) {




    /* Statistics */
    double startTime = 0.0, endTime = 0.0, avg, dev; /* Timing */
    double times[10]; /* Times for all runs */



    /* Setup MPI */


    int matrixDimensions;
    int nodes;

    initMPI(&argc, &argv);

    /* Get parameters */
    if (argc == 3) {
        /* Get number of processes */
        nodes = atoi(argv[1]);

        /* Get maximum matrix dimension */
        matrixDimensions = atoi(argv[2]);
    } else {
        printf("Wrong number of parameters\n");
        exit(-1);
    }

    /* Write header */
    if (rank == 0) {
        results_header();
    }

    /* Make cartesian grid */
    /*
    ...
    */

    /* Make and allocate matrices */


    double *matrixA = NULL;
    double *matrixB = NULL;
    double *matrixC = NULL;
//    allocate_matrices(matrixDimensions, matrixA,matrixB, matrixC);
    matrixA = malloc(matrixDimensions*matrixDimensions* sizeof(double));
    matrixB = malloc(matrixDimensions*matrixDimensions* sizeof(double));
    matrixC = malloc(matrixDimensions*matrixDimensions* sizeof(double));
    make_full_matrices(matrixDimensions,matrixA,matrixB);
    setzero(matrixDimensions, matrixC);


    /* Run each config 10 times */
    for (int k = 0; k < 10; k++) {

        /* Start timer */
        MPI_Barrier(MPI_COMM_WORLD);
        if (myrank == 0) {
            startTime = MPI_Wtime();
        }

        spread_matrix(nodes, matrixDimensions, matrixA, matrixB);
        if (coordinates[2] == 0) {
            BlockMatrixC = (double *) malloc(sizeof(double) * blockLength * blockLength);
        }



        /* Do work */
        distribute();
        broadcast();
        multiplyAndReduce();

        /* End timer */
        MPI_Barrier(MPI_COMM_WORLD);
        if (myrank == 0) {
            endTime = MPI_Wtime();
            times[k] = endTime - startTime;
        }
        /* Reset matrices */
        reset_matrices(matrixDimensions, matrixA, matrixB, matrixC);
    }

    /* Destroy matrices */
    destroy_matrices(matrixA,matrixB,matrixC);
    /* Print stats */


    if (myrank == 0) {
        dev = 0;
        for (int k = 0; k < 10; k++) {
//            printf("%f \n",times[k]);
        }
        avg = average(10, times, &dev);
        double Ts = (pow(matrixDimensions, 3)*2) / 16000000000;
//        printf("TS: %f \n",Ts);
        double Tp = avg;
//        printf("TP: %f \n",Tp);
        double S = Ts/Tp;
//        printf("S: %f \n",S);
        double eff = S / nodes;
        avg = average(10, times, &dev);
        write_result(matrixDimensions,nodes, avg, dev, eff);
    }

    /* Exit program */
    MPI_Finalize();

    return 0;
}


