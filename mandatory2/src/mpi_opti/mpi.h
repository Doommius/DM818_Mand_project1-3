//
// Created by jervelund on 11/18/17.
//
#include "common.h"
#ifndef DM818_MPI_H
#define DM818_MPI_H

#endif //DM818_MPI_H

typedef struct {
    int particleCount;
    particle_t *particles;
    int coordinateStart;
} edgezone;

typedef struct {
    int id;
    int step;
    double x;
    double y;
} localParticle;