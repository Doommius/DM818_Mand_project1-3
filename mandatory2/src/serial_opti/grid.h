//
// Created by jervelund on 11/12/17.
//

#include <vector>
#include "common.h"
#ifndef DM818_GRID_H
#define DM818_GRID_H



extern int grid_Columns;

void grid_init(int size);

void grid_add(particle_t particle);

void grid_move(particle_t particle);




#endif //DM818_GRID_H
