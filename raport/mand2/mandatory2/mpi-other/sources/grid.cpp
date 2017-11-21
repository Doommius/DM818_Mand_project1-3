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

void grid_clear_area(int start, int end){
    if (start >= 0 && end < grid.size()) {
        for (int i = start; i < end; i++) {
            grid[i].clear();
        }
    }
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
    int p_coordinate = grid_particle_index(particle);
//    printf("checking if cell is empty %s \n",grid[p_coordinate].empty() ? "true" : "false");
    grid[p_coordinate].clear();
//    printf("checking if cell is empty %s \n",grid[p_coordinate].empty() ? "true" : "false");

}

std::vector<particle_t *> grid_get(int n){
    return grid[n];
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

