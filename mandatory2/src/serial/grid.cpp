#include <vector>
#include <algorithm>
#include <assert.h>
#include <unordered_map>

#include "grid.h"

std::vector<std::vector<particle_t *> > grid;

int gridsize;

void grid_init(double size) {
    gridsize = (int) (size + 1);
    grid.resize((unsigned long) (gridsize * gridsize));
    printf("Grid size: %d\n", (int) grid.size());
}

int grid_GetParticleCoordinate(const particle_t *particle) {
    int coordinate = (int) (particle->x / 0.01) * gridsize + (int) (particle->y / 0.01);
    return coordinate;
}

void grid_Add(particle_t *particle) {
    // Idea is to see the array as the grid, and store particles into it
    // at their positions, converting to integers is nessecary (right?), and
    // should only slightly affect their real position so a few too many
    // might be included.

    // I fear the whole double positioning might have eluded me though
    // and this conversion wont work. Needs to be tested i guess.
    int coordinate = grid_GetParticleCoordinate(particle);
    grid[coordinate].push_back(particle);
}

void gridReset() {
    grid.clear();
    grid.resize((unsigned long) (gridsize * gridsize));
}

void gridValidate(int particle_count, const char *const file, int line) {
    int particles_in_system = 0;
    for (int i = 0; i < grid.size(); i++) {
        for (auto particle : grid[i]) {
            int realCoordinate = grid_GetParticleCoordinate(particle);
            if (realCoordinate != i) {
                printf("Assertion fail at %s:%d\n", file, line);
                printf("I found a particle in %d which belongs to %d\n", i, realCoordinate);
                assert(realCoordinate == i);
            }
        }
        particles_in_system += grid[i].size();
    }
    if (particle_count != 0) {
        if (particles_in_system != particle_count) {
            printf("Assertion fail at %s:%d\n", file, line);
            assert(particles_in_system == particle_count);
        }
    }
}

void gridValidateParticlesOnlyWithinZone(int startInclusive, int endInclusive) {
    assert(startInclusive >= 0);
    assert(endInclusive <= gridsize * gridsize);

//    for (int i = 0; i < startInclusive; i++) {
//        CASSERT(grid[i].empty(), "Expected cell %d to be empty. We're only expecting particles in %d to %d", i,
//                startInclusive, endInclusive);
//    }
//    for (int i = endInclusive + 1; i < gridColumns * gridColumns; i++) {
//        CASSERT(grid[i].empty(), "Expected cell %d to be empty. We're only expecting particles in %d to %d", i,
//                startInclusive, endInclusive);
//    }
}


void grid_Remove(particle_t *particle) {
    int coordinate = grid_GetParticleCoordinate(particle);

    std::vector<particle_t *> &cell = grid[coordinate];
    cell.erase(std::remove(cell.begin(), cell.end(), particle), cell.end());
}

void gridPurge(int startInclusive, int endExclusive) {
    if (startInclusive >= 0 && endExclusive < grid.size()) {
        for (int i = startInclusive; i < endExclusive; i++) {
            grid[i].clear();
        }
    }
}

std::vector<particle_t *> gridGetAt(int coordinate) {
    return grid[coordinate];
}

std::vector<particle_t *> gridGetCollisionsAtNeighbor(particle_t *particle, int offsetX, int offsetY) {
    int coordinate = grid_GetParticleCoordinate(particle);
    int x = (coordinate % gridsize) + offsetX;
    int y = (coordinate / gridsize) + offsetY;

    if (x < 0 || y < 0 || x >= gridsize || y >= gridsize) {
        return std::vector<particle_t *>();
    }

    return grid[coordinate + (offsetY * gridsize) + offsetX];
}
