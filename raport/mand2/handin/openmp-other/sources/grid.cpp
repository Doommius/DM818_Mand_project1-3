//
// Created by jervelund on 11/12/17.
//

#include <cstdio>
#include <algorithm>
#include <vector>
#include <cmath>
#include "grid.h"

std::vector<std::vector<particle_t *> > grid;

unsigned int grid_size;


void grid_init(int n) {
    grid_size = (unsigned int) (sqrt(0.0005 * n) / 0.01) + 1;
    grid.resize((grid_size * grid_size));
}

/*
 * useless
 */
int grid_get_size() {
    int returnval = grid_size;
    return returnval;
}


int grid_particle_index(particle_t *particle) {

    return (int) (particle->x / 0.01) * grid_size + (int) (particle->y / 0.01);
}

/**
 * we see the grid as a matrix and insert it, in this way we can easily split the array up into blocks.
 * @param particle
 */
void grid_add(particle_t *particle) {
    grid[grid_particle_index(particle)].push_back(particle);
}

/**
 * Clears the whole grid.
 *
 */
void grid_purge() {
    grid.clear();
    grid.resize(grid_size * grid_size);
}

void grid_remove(particle_t *particle) {
    int coordinate = grid_particle_index(particle);
    std::vector<particle_t *> &cell = grid[coordinate];
    cell.erase(std::remove(cell.begin(), cell.end(), particle), cell.end());
}
std::vector<particle_t *> gridGetCollisionsAtNeighbor(particle_t *particle, int offsetX, int offsetY) {
    int p_cord = grid_particle_index(particle);
    int x = (p_cord % grid_size) + offsetX;
    int y = (p_cord / grid_size) + offsetY;
    if (x < 0 || y < 0 || x >= grid_size || y >= grid_size) {
        return std::vector<particle_t *>();
    }
    return grid[p_cord + (offsetY * grid_size) + offsetX];
}

