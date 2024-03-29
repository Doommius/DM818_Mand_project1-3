//
// Created by jervelund on 11/11/17.
//
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <algorithm>
#include <assert.h>
#include <unordered_map>
#include <ostream>
#include <iostream>

#include "common.h"
#include "grid.h"
#include "openmp_opti.h"

// Used to account for time spent
std::unordered_map<std::string, double> profiling;

//
//  benchmarking program
//
int maxRows = 0;

int main(int argc, char **argv) {
    int navg, nabsavg=0;
    double dmin, absmin=1.0,davg,absavg=0.0;
    double rdavg,rdmin;
    int rnavg;

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

    n = 1000;

    char *savename = read_string(argc, argv, "-o", NULL);
    char *sumname = read_string(argc, argv, "-s", NULL);

    FILE *fsave = savename ? fopen(savename, "w") : NULL;
    FILE *fsum = sumname ? fopen(sumname, "a") : NULL;

    particle_t *particles = (particle_t *) malloc(n * sizeof(particle_t));
    set_size(n);
    init_particles(n, particles);


    /**
     * setup and add elements to grid, move move to its own function.
     */
    grid_init(sqrt(0.0005 * n) / 0.01);

    for (int i = 0; i < n; i++) {
        grid_Add(&particles[i]);
    }
    //
    //  simulate a number of time steps
    //


    double simulation_time = read_timer();
    for (int step = 0; step < NSTEPS; step++) {
        //
        //  compute forces
        //
        for (int i = 0; i < n; i++) {
            particles[i].ax = particles[i].ay = 0;

            for (int row = -1; row < 1; row++) {
                for (int column = -1; column <= 1; column++) {
                    const std::vector<particle_t *> &cell =
                            gridGetCollisionsAtNeighbor(&particles[i], row, column);
                    for (auto particle : cell) {
                        apply_force( particles[i], *particle, &dmin, &davg, &navg );
                    }

                }
            }

        }


        if( find_option( argc, argv, "-no" ) == -1 )
        {
            //
            // Computing statistical data
            //
            if (navg) {
                absavg +=  davg/navg;
                nabsavg++;
            }
            if (dmin < absmin) absmin = dmin;

            //
            //  save if necessary
            //
            if( fsave && (step%SAVEFREQ) == 0 )
                save( fsave, n, particles );
        }
        //
        //  move particles
        //
        for (int i = 0; i < n; i++) {
            int coordinate = grid_GetParticleCoordinate(&particles[i]);
            grid_Remove(&particles[i]);
            move(particles[i]);
            grid_Add(&particles[i]);
            int newCoordinate = grid_GetParticleCoordinate(&particles[i]);
            if (abs(coordinate - newCoordinate) > gridsize) {
                int i1 = abs(coordinate - newCoordinate) / gridsize;
//                printf("Moved %d cells (From %d to %d [%d rows])\n", abs(coordinate - newCoordinate), coordinate,
//                       newCoordinate, i1);
                if (i1 > maxRows) maxRows = i1;
            }
        }

        //  save if necessary
        if (fsave && (step % SAVEFREQ) == 0) {
//            BEGIN_TIMED_ZONE(saveLocalResult);
            save(fsave, n, particles);
//            END_TIMED_ZONE(saveLocalResult);
        }
    }
    simulation_time = read_timer() - simulation_time;

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
        fprintf(fsum, "%d %g\n", n, simulation_time);

    //
    // Clearing space
    //
    if (fsum)
        fclose(fsum);
    free(particles);
    if (fsave)
        fclose(fsave);

    return 0;
}
