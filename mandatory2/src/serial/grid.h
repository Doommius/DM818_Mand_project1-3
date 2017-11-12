//
// Created by jervelund on 11/11/17.
//

#include <vector>
#include "common.h"

extern int gridsize;

int grid_GetParticleCoordinate(const particle_t *particle);

void grid_init(double size);

void grid_Add(particle_t *particle);

void grid_Remove(particle_t *particle);

void grid_Validate(int particle_count, const char *const file, int line);

std::vector<particle_t *> gridGetCollisionsAtNeighbor(particle_t *particle, int offsetX, int offsetY);



std::vector<particle_t *> gridGetAt(int coordinate);

void grid_Reset();

void grid_Purge(int startInclusive, int endExclusive);

void grid_ValidateParticlesOnlyWithinZone(int startInclusive, int endInclusive);
