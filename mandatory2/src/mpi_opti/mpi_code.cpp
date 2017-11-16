//
// Created by jervelund on 11/16/17.
//

#include <cstring>
#include <mpi.h>
#include "mpi_code.h"
#include "grid.h"


MPI_Datatype particleStructure;

void prepareEdge(edgezone &zone) {
    // Move the particles from the real grid into the buffer
    zone.particleCount = 0;
    for (int i = zone.coordinateStart; i < zone.coordinateStart + grid_get_size(); i++) {
        for (auto particle : grid_get(i)) {
            memcpy(&zone.particles[zone.particleCount], particle, sizeof(particle_t));
            zone.particleCount++;
        }
    }
}

void exchangeEdge(edgezone &local, edgezone &remote, int multiplier,int rank, MPI_Datatype particleType) {
    // Communicate number of particles in update
    MPI_Sendrecv(&local.particleCount, 1, MPI_INT, rank + (1 * multiplier), 0, &remote.particleCount, 1,
                 MPI_INT, rank + (1 * multiplier), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    // Exchange particles
    MPI_Sendrecv(local.particles, local.particleCount, particleType, rank + (1 * multiplier), 0,
                 remote.particles, remote.particleCount, particleType, rank + (1 * multiplier), 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);

}


void initParticleType() {
    int count = 7;
    int blocklens[] = {1, 1, 1, 1, 1, 1, 1};

    MPI_Aint indices[7];
    indices[0] = (MPI_Aint) offsetof(particle_t, x);
    indices[1] = (MPI_Aint) offsetof(particle_t, y);
    indices[2] = (MPI_Aint) offsetof(particle_t, vx);
    indices[3] = (MPI_Aint) offsetof(particle_t, vy);
    indices[4] = (MPI_Aint) offsetof(particle_t, ax);
    indices[5] = (MPI_Aint) offsetof(particle_t, ay);

    MPI_Datatype old_types[] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_INT};

    MPI_Type_struct(count, blocklens, indices, old_types, &particleStructure);
    MPI_Type_commit(&particleStructure);
}