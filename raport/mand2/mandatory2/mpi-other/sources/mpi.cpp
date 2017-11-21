#include <mpi.h>
#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <cmath>
#include "common.h"
#include "grid.h"
#include "mpi.h"

//
//  benchmarking program
//

int cellsperthread;
int rank;
int maxrank;
int totalParticleCount;
int localparticlescount;

edgezone localupper;
edgezone locallower;
edgezone remoteupper;
edgezone remotelower;

MPI_Datatype PARTICLE;

std::vector<particle_t> insertionsIntoUpperBorrowed;
std::vector<particle_t> insertionsIntoLowerBorrowed;

//int maxPosition;
particle_t *localParticles;


void prepareEdge(edgezone &edge) {
    // Move the particles from the real grid into the buffer
    edge.particleCount = 0;
    for (int i = edge.coordinateStart; i < edge.coordinateStart + grid_get_size(); i++) {
        for (auto particle : grid_get(i)) {
            memcpy(&edge.particles[edge.particleCount], particle, sizeof(particle_t));
            edge.particleCount++;
        }
    }
}

void exchangeEdge(edgezone &local, edgezone &remote, int multiplier) {
    // Communicate number of particles in update
    MPI_Sendrecv(&local.particleCount, 1, MPI_INT, rank + (1 * multiplier), 0, &remote.particleCount, 1,
                 MPI_INT, rank + (1 * multiplier), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    // Exchange particles
    MPI_Sendrecv(local.particles, local.particleCount, PARTICLE, rank + (1 * multiplier), 0,
                 remote.particles, remote.particleCount, PARTICLE, rank + (1 * multiplier), 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);

}


particle_t *exchangeInfomation(std::vector<particle_t> &particlesToExchange, int multiplier, int *recievedCount) {
    int count;
    auto sendingCount = (int) particlesToExchange.size();


    MPI_Sendrecv(&sendingCount, 1, MPI_INT, rank + (1 * multiplier), 0, &count, 1, MPI_INT, rank + (1 * multiplier), 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    int i = 0;
    auto *prepared = (particle_t *) malloc(sizeof(particle_t) * particlesToExchange.size());
    for (auto particle : particlesToExchange) {
        memcpy(&prepared[i], &particle, sizeof(particle_t));
        i++;
    }


    particle_t *receivedParticles;
    receivedParticles = (particle_t *) malloc(sizeof(particle_t) * count);

    MPI_Sendrecv(prepared, sendingCount, PARTICLE, rank + (1 * multiplier), 0, receivedParticles, count, PARTICLE,
                 rank + (1 * multiplier), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    *recievedCount = count;
    free(prepared);
    return receivedParticles;
}


void exchangeInformationAbove(particle_t **insertedLower, int *insertedLowerCount) {
    if (rank < maxrank - 1) {
        exchangeEdge(localupper, remoteupper, 1);
        *insertedLower = exchangeInfomation(insertionsIntoUpperBorrowed, 1, insertedLowerCount);

    }
}

void exchangeInformationBelow(particle_t **insertedUpper, int *insertedUpperCount) {
    if (rank > 0) {
        exchangeEdge(locallower, remotelower, -1);
        *insertedUpper = exchangeInfomation(insertionsIntoLowerBorrowed, -1, insertedUpperCount);
    }
}

void exchangeInformation(){
    int insertedIntoUpperOwnedCount = 0;
    particle_t *insertedIntoUpperOwned = nullptr;
    int insertedIntoLowerOwnedCount = 0;
    particle_t *insertedIntoLowerOwned = nullptr;

    prepareEdge(localupper);
    prepareEdge(locallower);

    grid_clear_area(remotelower.coordinateStart,remotelower.coordinateStart+grid_get_size());
    grid_clear_area(remoteupper.coordinateStart,remoteupper.coordinateStart+grid_get_size());

    if (rank % 2 == 0) {
        exchangeInformationAbove(&insertedIntoLowerOwned, &insertedIntoLowerOwnedCount);
        exchangeInformationBelow(&insertedIntoUpperOwned, &insertedIntoUpperOwnedCount);
    }else {
        exchangeInformationBelow(&insertedIntoUpperOwned, &insertedIntoUpperOwnedCount);
        exchangeInformationAbove(&insertedIntoLowerOwned, &insertedIntoLowerOwnedCount);
    }


}

int main(int argc, char **argv) {
    int navg, nabsavg = 0;
    double dmin, absmin = 1.0, davg, absavg = 0.0;
    double rdavg, rdmin;
    int rnavg;

    //
    //  process command line parameters
    //
    if (find_option(argc, argv, "-h") >= 0) {
        printf("Options:\n");
        printf("-h to see this help\n");
        printf("-n <int> to set the number of particles\n");
        printf("-o <filename> to specify the output file name\n");
        printf("-s <filename> to specify a summary file name\n");
        printf("-no turns off all correctness checks and particle output\n");
        return 0;
    }

    int n = read_int(argc, argv, "-n", 1000);
    char *savename = read_string(argc, argv, "-o", nullptr);
    char *sumname = read_string(argc, argv, "-s", nullptr);

    //
    //  set up MPI
    //
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &maxrank);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);


    MPI_Type_contiguous(6, MPI_DOUBLE, &PARTICLE);
    MPI_Type_commit(&PARTICLE);





    //
    //  allocate generic resources
    //
    FILE *fsave = savename && rank == 0 ? fopen(savename, "w") : nullptr;
    FILE *fsum = sumname && rank == 0 ? fopen(sumname, "a") : nullptr;


    auto *particles = (particle_t *) malloc(n * sizeof(particle_t));
    localParticles = (particle_t *) malloc(totalParticleCount * sizeof(particle_t));

    auto *particlesToSend = (particle_t *) malloc(totalParticleCount * sizeof(particle_t));

    int sendCount[maxrank];
    int sendDisplacement[maxrank];


    grid_init(totalParticleCount);
    set_size(totalParticleCount);

    cellsperthread = grid_get_size() * (int) ceil(grid_get_size() / (double) maxrank);


    if (rank == 0) {
        int sendStart = 0;
        int currentProcessor = 0;
        int counter = 0;

        init_particles(totalParticleCount, particles);

        for (int i = 0; i < totalParticleCount; ++i) {
            grid_add(&particles[i]);
        }
        for (int i = 0; i < grid_get_size() * grid_get_size(); i++) {
            auto cell = grid_get(i);
            for (auto particle : cell) {
                memcpy(particlesToSend + counter, particle, sizeof(particle_t));
                counter++;
            }
            if ((i > 0 && i % cellsperthread == 0) || i == (grid_get_size() * grid_get_size()) - 1) {
                sendCount[currentProcessor] = counter - sendStart;
                sendDisplacement[currentProcessor] = sendStart;
                sendStart = counter;
                currentProcessor++;
            }
        }
        grid_purge();

    }

    // Distribute the particle count to all processors
    MPI_Scatter(sendCount, 1, MPI_INT, localParticles, 1, MPI_INT, 0, MPI_COMM_WORLD);
    // Distribute the actual particles
    MPI_Scatterv(particlesToSend, sendCount, sendDisplacement, PARTICLE, localParticles, totalParticleCount,
                 PARTICLE,
                 0, MPI_COMM_WORLD);

    remotelower.particles = (particle_t *) malloc(sizeof(particle_t) * totalParticleCount);
    remoteupper.particles = (particle_t *) malloc(sizeof(particle_t) * totalParticleCount);
    remotelower.coordinateStart = (cellsperthread * rank) - grid_get_size();
    remoteupper.coordinateStart = (cellsperthread * (rank + 1));

    locallower.particles = (particle_t *) malloc(sizeof(particle_t) * totalParticleCount);
    localupper.particles = (particle_t *) malloc(sizeof(particle_t) * totalParticleCount);
    locallower.coordinateStart = cellsperthread * rank;
    localupper.coordinateStart = min((cellsperthread * (rank + 1)) - grid_get_size(), grid_get_size() * (grid_get_size() - 1));


    for (int i = 0; i < localparticlescount; i++) {
        grid_add(&localParticles[i]);
    }

    free(particles);
    free(particlesToSend);

    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer();
    for (int step = 0; step < NSTEPS; step++) {
        navg = 0;
        dmin = 1.0;
        davg = 0.0;
        //
        //  save current step if necessary (slightly different semantics than in other codes)
        //
        if (find_option(argc, argv, "-no") == -1)
            if (fsave && (step % SAVEFREQ) == 0)
                save(fsave, n, particles);

        MPI_Barrier(MPI_COMM_WORLD);
        exchangeInformation();


        //
        //  compute all forces
        //
        for (int i = 0; i < localparticlescount; i++) {
            // traverse included neighbors
            for (int offsetX = -1; offsetX <= 1; offsetX++) {
                for (int offsetY = -1; offsetY <= 1; offsetY++) {
                    const std::vector<particle_t *> &cell =
                            gridGetCollisionsAtNeighbor(&localParticles[i], offsetX, offsetY);
                    for (auto particle : cell) {
                        apply_force(localParticles[i], *particle, &dmin, &davg, &navg);
                    }
                }
            }
        }

        if (find_option(argc, argv, "-no") == -1) {

            MPI_Reduce(&davg, &rdavg, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            MPI_Reduce(&navg, &rnavg, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
            MPI_Reduce(&dmin, &rdmin, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);


            if (rank == 0) {
                //
                // Computing statistical data
                //
                if (rnavg) {
                    absavg += rdavg / rnavg;
                    nabsavg++;
                }
                if (rdmin < absmin) absmin = rdmin;
            }
        }

        //
        //  move particles
        //
        grid_purge();
        for (int i = 0; i < localparticlescount; i++){

            move(localParticles[i]);
            int cord = grid_particle_index(&localParticles[i]);
            if (cord < remotelower.coordinateStart + grid_get_size() ||
                    cord >= remoteupper.coordinateStart) { // Out of bounds
                // Copy the particle and insert into a vector of insertions into a certain ghost zone
                particle_t copy;
                memcpy(&copy, &localParticles[i], sizeof(particle_t));

                if (cord >= remotelower.coordinateStart &&
                        cord < remotelower.coordinateStart + grid_get_size()) {
                    insertionsIntoLowerBorrowed.push_back(copy);
                } else if (cord >= remoteupper.coordinateStart &&
                        cord < remoteupper.coordinateStart + grid_get_size()) {
                    insertionsIntoUpperBorrowed.push_back(copy);
                }

            } else {
                grid_add(&localParticles[i]);
            }
        }


    }
    simulation_time = read_timer() - simulation_time;

    if (rank == 0) {
        printf("n = %d, simulation time = %g seconds", n, simulation_time);

        if (find_option(argc, argv, "-no") == -1) {
            if (nabsavg) absavg /= nabsavg;
            //
            //  -The minimum distance absmin between 2 particles during the run of the simulation
            //  -A Correct simulation will have particles stay at greater than 0.4 (of cutoff) with typical values between .7-.8
            //  -A simulation where particles don't interact correctly will be less than 0.4 (of cutoff) with typical values between .01-.05
            //
            //  -The average distance absavg is ~.95 when most particles are interacting correctly and ~.66 when no particles are interacting
            //
            printf(", absmin = %lf, absavg = %lf", absmin, absavg);
            if (absmin < 0.4) printf("\nThe minimum distance is below 0.4 meaning that some particle is not interacting");
            if (absavg < 0.8) printf("\nThe average distance is below 0.8 meaning that most particles are not interacting");
        }
        printf("\n");

        //
        // Printing summary data
        //
        if (fsum)
            fprintf(fsum, "%d %d %g\n", n, maxrank, simulation_time);
    }

    //
    //  release resources
    //
    if (fsum)
        fclose(fsum);
    free(localParticles);
    free(particles);
    if (fsave)
        fclose(fsave);

    MPI_Finalize();

    return 0;
}
