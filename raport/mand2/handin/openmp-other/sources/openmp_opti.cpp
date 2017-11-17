//#include <stdlib.h>
//#include <stdio.h>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <unordered_map>

#include "common.h"
#include "grid.h"
#include "omp.h"

//double starttimer;
//double endtimer;
//double inittimer = 0;
//double movetimer = 0;
//double applyforcetimer = 0;

//
//  benchmarking program
//
int main(int argc, char **argv) {

    int navg, nabsavg = 0, numthreads = 0;
    double davg, dmin, absmin = 1.0, absavg = 0.0;

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
//    n = 100;

    char *savename = read_string(argc, argv, "-o", NULL);
    char *sumname = read_string(argc, argv, "-s", NULL);

    FILE *fsave = savename ? fopen(savename, "w") : NULL;
    FILE *fsum = sumname ? fopen(sumname, "a") : NULL;

    auto *particles = (particle_t *) malloc(n * sizeof(particle_t));
    set_size(n);
    init_particles(n, particles);

    //
    // Setup and Add particles to grid
    //
//    starttimer = read_timer();
    grid_init(n);

    for (int i = 0; i < n; i++) {
        grid_add(&particles[i]);
    }
//    endtimer = read_timer();
//    inittimer = endtimer - starttimer;
    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer();

#pragma omp parallel
    {
        if (omp_get_thread_num() == 0) {
            numthreads = omp_get_num_threads();
        };

        for (int step = 0; step < NSTEPS; step++) {
//            printf("at step %i \n",step);

            navg = 0;
            davg = 0.0;
            dmin = 1.0;
            //
            //  compute forces
            //
//            starttimer = read_timer();
//            printf("apply force \n");
#pragma omp barrier
#pragma omp for
            for (int i = 0; i < n; i++) {
                particles[i].ax = particles[i].ay = 0;
                for (int j = -1; j <= 1; j++) {
                    for (int k = -1; k <= 1; k++) {
                        const std::vector<particle_t *> &cell = gridGetCollisionsAtNeighbor(&particles[i], j, k);
                        for (particle_t *particle : cell) {
                            apply_force(particles[i], *particle, &dmin, &davg, &navg);

                        }


                    }

                }

            }
//            endtimer = read_timer();
//            applyforcetimer += endtimer - starttimer;

            //
            //  move particles
            // Seems to work correctly when using purge, but breaks with using remove.
            //

//            starttimer = read_timer();


            if (omp_get_thread_num() == 0) {
                grid_purge();
            }
#pragma omp barrier
            for (int i = 0; i < n; i++) {
//                grid_remove(&particles[i]);

                move(particles[i]);
                grid_add(&particles[i]);

            }
//            endtimer = read_timer();
//            movetimer += endtimer - starttimer;


            if (find_option(argc, argv, "-no") == -1) {
                //
                // Computing statistical data
                //
#pragma omp master
                if (navg) {
                    absavg += davg / navg;
                    nabsavg++;
                }
#pragma omp critical
                if (dmin < absmin) absmin = dmin;

                //
                //  save if necessary
                //
#pragma omp master
                if (fsave && (step % SAVEFREQ) == 0)
                    save(fsave, n, particles);
            }
        }
    }
    simulation_time = read_timer() - simulation_time;

    printf("n = %d,threads = %d, simulation time = %g seconds", n, numthreads, simulation_time);
//    printf("Init timer = %f \n", inittimer);
//    printf("applyforce timer = %f \n", applyforcetimer);
//    printf("move timer = %f \n", movetimer);

    if (find_option(argc, argv, "-no") == -1) {
        if (nabsavg) absavg /= nabsavg;
        //
        //  -The minimum distance absmin between 2 particles during the run of the simulation
        //  -A Correct simulation will have particles stay at greater than 0.4 (of cutoff) with typical values between .7-.8
        //  -A simulation where particles don't interact correctly will be less than 0.4 (of cutoff) with typical values between .01-.05
        //
        //  -The average distance absavg is ~.95 when most particles are interacting correctly and ~.66 when no particles are interacting
        //
        printf("absmin = %lf, absavg = %lf", absmin, absavg);
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
