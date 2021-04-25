#pragma once
#include <stdio.h>
#include <vector>
#define STB_IMAGE_IMPLEMENTATION
#include "stb_read_image.h"

// TODO: please better names please please please
struct {
    double **pot; // potential: Device and Host
    double *density; // charge density: Host only
    int *conductor; // condutors with id: Device and host
    unsigned char ci; // current grid index
    int width, height;
    
    // Information about the conductors
    double *charges; // Device only
    double *capacity; // Device only TODO: already use the matrix
    int n_conductors;
} typedef Grid;

int in_bounds(Grid *grid, int x, int y) {
    return (x < grid->width && x >= 0) && (y < grid->height && y >= 0);
}

int get_index(Grid *grid, int x, int y) {
    return y * grid->width + x;
}

void init_grid(Grid* grid, int width, int height) {
    grid->width = width;
    grid->height = height;
    grid->ci = 0;

    cudaMallocManaged(&grid->pot, 2 * sizeof(double*));
    
    cudaMallocManaged(&grid->conductor, width * height * sizeof(int));
    cudaMallocManaged(&grid->pot[0], width * height * sizeof(double));
    cudaMallocManaged(&grid->pot[1], width * height * sizeof(double));
    cudaMallocManaged(&grid->density, width * height * sizeof(double));
}

void copy_zoomed(Grid* grid, Grid* ref, int zoom) {
    init_grid(grid, ref->width*zoom, ref->height*zoom);
    
    // Copy all the information zoomed
    for (int x = 0; x < grid->width; x++) {
        for (int y = 0; y < grid->height; y++) {
            int g_index = get_index(grid, x, y);
            int r_index = get_index(ref, x/zoom, y/zoom);

            grid->pot[0][g_index] = ref->pot[0][r_index];
            grid->pot[1][g_index] = ref->pot[1][r_index];

            grid->conductor[g_index] = ref->conductor[r_index];
            grid->density[g_index] = ref->density[r_index];
        }
    }

    // Initialize all the arrays
    cudaMallocManaged(&grid->capacity, ref->n_conductors * ref->n_conductors * sizeof(double));
    cudaMallocManaged(&grid->charges, ref->n_conductors * sizeof(double));
    grid->n_conductors = ref->n_conductors;

    for (int i = 0; i < ref->n_conductors; i++) {
        grid->charges[i] = ref->charges[i];
    }

    grid->ci = ref->ci;
}

// Recursive fill algorithm
void fill_with(Grid* grid, unsigned char* img, int x, int y, int color) {
    int index = get_index(grid, x, y);
    grid->conductor[index] = color;

    int dirs[4][2] = { {-1, 0}, {0, -1}, {1, 0}, {0, 1} };
    for (int i = 0; i < 4; i++) {
        int nx = x + dirs[i][0];
        int ny = y + dirs[i][1];

        if (!in_bounds(grid, nx, ny)) {
            continue;
        }

        int n_index = get_index(grid, nx, ny);
        if (grid->conductor[n_index] == 0 && img[n_index] == img[index]) {
            fill_with(grid, img, nx, ny, color);
        }
    }
}


// Neutral charge is 127
void from_image(Grid* grid, unsigned char* img, int width, int height) {
    int c_count = 0;

    init_grid(grid, width, height);

    std::vector<int> cond_charges;
    for (int x = 0; x < width; x++) {
        for (int y = 0; y < height; y++) {
            int index = get_index(grid, x, y);
            if (img[index] != 0 && grid->conductor[index] == 0) {
                fill_with(grid, img, x, y, c_count+1);
                cond_charges.push_back((int)img[index] - 127);
                
                c_count++;
            }
        }
    }

    // Initialize all the arrays
    cudaMallocManaged(&grid->capacity, c_count * c_count * sizeof(double));
    cudaMallocManaged(&grid->charges, c_count * sizeof(double));
    grid->n_conductors = c_count;

    // Copy from vector to pointer
    for (int i = 0; i < c_count; i++) {
        grid->charges[i] = cond_charges[i];
        printf("# %d %d\n", i, cond_charges[i]);
    }
}

