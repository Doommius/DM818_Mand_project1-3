//
// Created by jervelund on 11/12/17.
//

#include <cstdio>
#include <vector>
#include "grid.h"

extern int grid_Columns;

std::vector<std::vector<particle_t *> > grid;


void grid_init(double size){
    grid_Columns = (int) (size + 1);
    grid.resize((unsigned long) (grid_Columns * grid_Columns));
    printf("Grid size: %d\n", (int) grid.size());
}


int grid_particle_index(particle_t *particle){

    return (int) (particle->x / 0.01) * grid_Columns + (int) (particle->y / 0.01);
}

/**
 * we see the grid as a matrix and insert it, in this way we can easily split the array up into blocks.
 * @param particle
 */
void grid_Add(particle_t *particle) {
    grid[grid_particle_index(particle)].push_back(particle);
}

void grid_Move(particle_t *particle){
    std::vector<particle_t *> &cell = grid[grid_particle_index(particle)];
    cell.erase(std::remove(cell.begin(), cell.end(), particle), cell.end());
    move(reinterpret_cast<particle_t &>(particle));
    grid[grid_particle_index(particle)].push_back(particle);
}

particle_t * gridGetAt(int coordinate);

particle_t * gridGetAt(int coordinate) {
    return reinterpret_cast<particle_t *>(grid[coordinate].data());
}
