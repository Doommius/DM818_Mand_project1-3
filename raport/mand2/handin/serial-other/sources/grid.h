//
// Created by jervelund on 11/12/17.
//

#include <vector>
#include "common.h"
#ifndef DM818_GRID_H
#define DM818_GRID_H


int grid_get_size();

void grid_init(int size);

void grid_clear();

void grid_add(particle_t *particle);

void grid_remove(particle_t *particle);

std::vector<particle_t *> gridGetNeighbors(particle_t *particle, int indexx, int indexy);



#endif //DM818_GRID_H
