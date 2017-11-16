//
// Created by jervelund on 11/12/17.
//

#include <vector>
#include "common.h"
#ifndef DM818_GRID_H
#define DM818_GRID_H


int grid_get_size();

void grid_init(int size);

void grid_purge();

std::vector<particle_t *> grid_get(int n);

void grid_add(particle_t *particle);

void grid_remove(particle_t *particle);

std::vector<particle_t *> gridGetCollisionsAtNeighbor(particle_t *particle, int offsetX, int offsetY);



#endif //DM818_GRID_H
