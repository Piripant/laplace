#pragma once
#include <stdio.h>
#include "grid.cu"
#define STB_IMAGE_IMPLEMENTATION

__global__ void laplace(Grid *grid) {
    int x = threadIdx.x;
    int y = blockIdx.x;

    //int index = y * blockDim.x + x;
    int index = y * grid->width + x;

    // Only update the free space around the conductors
    if (grid->conductor[index] == 0) {
        double v = 0;
        
        if (y+1 < grid->height) {
            v += grid->pot[grid->ci][(y+1) * grid->width + x];
        }
        if (y-1 >= 0) {
            v += grid->pot[grid->ci][(y-1) * grid->width + x];
        }
        if (x+1 < grid->width) {
            v += grid->pot[grid->ci][y * grid->width + x+1];
        }
        if (x-1 >= 0) {
            v += grid->pot[grid->ci][y * grid->width + x-1];
        }

        // All the non explored potentials count as 0
        grid->pot[!grid->ci][index] = v / 4.0;
    }
}

// Compute the potential using an parrallel Jacobi method with laplace equation
int compute_potential(Grid *grid, double threshold) {    
    int steps = 0;
    double diff = threshold;
    while (diff >= threshold) {
        laplace<<<grid->width, grid->height>>>(grid);
        cudaDeviceSynchronize();
        grid->ci = !grid->ci; // The future is now the present

        // Compute the difference only every 1024 iterations
        if (steps % 1024 == 0) {
            diff = 0.0;
            // Compute the biggest difference between the two grids
            for (int i = 0; i < grid->width * grid->height; i++) {
                double delta = fabs(grid->pot[grid->ci][i] - grid->pot[!grid->ci][i]);
                
                if (delta > diff) {
                    diff = delta;
                }
            }
        }

        steps++;
    }

    return steps;
}

// Compute the density using Columbs Theorem
void compute_density(Grid* grid) {    
    int dirs[4][2] = { {0, 1}, {0, -1}, {1, 0}, {-1, 0} };
    for (int x = 1; x < grid->width - 1; x++) {
        for (int y = 1; y < grid->height - 1; y++) {
            int index = get_index(grid, x, y);

            if (grid->conductor[index] != 0) {
                for (int i = 0; i < 4; i++) {
                    int n_index = get_index(grid, x + dirs[i][0], y + dirs[i][1]);
                    
                    // The other block is not inside the conductor
                    if (grid->conductor[n_index] != grid->conductor[index]) {
                        // Approximation of the partial derivative in that direction of the field
                        double field = (grid->pot[grid->ci][index] - grid->pot[grid->ci][n_index]);
                        grid->density[index] = field;
                        break;
                    }
                }
            }
        }
    }
}


#include <armadillo>
using namespace arma;

// TODO: Make this into a tidier version with maybe option to print on files and file names
// TODO: Experimento with full induction of spherical conductors one inside the other
// Compute the capacities of each conductor
void compute_capacities(Grid* grid, double threshold) {
    double default_pot = 1.0;

    for (int n = 0; n < grid->n_conductors; n++) {
        printf("Resetting %d\n", n);

        // Reset the grid and set all the parts of this conductor to potential pot
        for (int i = 0; i < grid->width * grid->height; i++) {
            grid->density[i] = 0.0;
            if (grid->conductor[i] == n+1) {
                grid->pot[grid->ci][i] = default_pot;
                grid->pot[!grid->ci][i] = default_pot;
            } else {
                grid->pot[grid->ci][i] = 0.0;
                grid->pot[!grid->ci][i] = 0.0;
            }
        }

        printf("Computing %d\n", n);

        // Compute charge density and total charge
        int steps = compute_potential(grid, threshold);
        printf("Converged after %d\n", steps);
        //printf("Computing density %d\n", n);
        
        compute_density(grid);
        // Fill a row of the capacities matrix
        for (int c = 0; c < grid->n_conductors; c++) {
            double charge = 0.0;
            for (int i = 0; i < grid->width * grid->height; i++) {
                if (grid->conductor[i] == c+1) {
                    charge += grid->density[i];
                }
            }

            //printf("potential %d againts %d: %lf\n", n, c, charge / default_pot);
            grid->capacity[n + c*grid->n_conductors] = charge / default_pot;
        }

        // Print the potential
        char file_name[20];
        sprintf(file_name, "test_cond/potential%d.dat", n);
        FILE* pot_fp = fopen(file_name, "w");
        if (pot_fp == NULL) {
            fprintf(stderr, "Error loading potential output file");
        }
        for (int x = 0; x < grid->width; x++) {
            for (int y = 0; y < grid->height; y++) {
                int index = get_index(grid, x, y);
                fprintf(pot_fp, "%d %d %lf\n", x, y, grid->pot[grid->ci][index]);
            }
        }
        
        // Print the density
        char den_name[20];
        sprintf(den_name, "test_cond/density%d.dat", n);
        FILE* den_fp = fopen(den_name, "w");
        if (den_fp == NULL) {
            fprintf(stderr, "Error loading density output file");
        }
        for (int x = 0; x < grid->width; x++) {
            for (int y = 0; y < grid->height; y++) {
                int index = get_index(grid, x, y);
                fprintf(den_fp, "%d %d %lf\n", x, y, grid->density[index]);
            }
        }
    }
    
    // Invert the matrix
    mat A(grid->n_conductors, grid->n_conductors);
    for (int i = 0; i < grid->n_conductors; i++) {
        for (int j = 0; j < grid->n_conductors; j++) {
            A(i, j) = grid->capacity[i + j*grid->n_conductors];
        }
    }

    A.print("Capacitance: ");

    mat B = inv(A);
    B.print("inv(A): ");

    // Assign potentials
    double* pots = (double*)calloc(grid->n_conductors, sizeof(double));
    for (int i = 0; i < grid->n_conductors; i++) {
        for (int j = 0; j < grid->n_conductors; j++) {
            pots[i] += grid->charges[j] * B(i, j);
        }
        printf("Potential of %d: %lf\n", i, pots[i]);
    }

    // Set grid
    for (int i = 0; i < grid->width * grid->height; i++) {
        grid->density[i] = 0.0;

        if (grid->conductor[i] != 0) {
            grid->pot[grid->ci][i] = pots[grid->conductor[i]-1];
            grid->pot[!grid->ci][i] = pots[grid->conductor[i]-1];
        } else {
            grid->pot[grid->ci][i] = 0.0;
            grid->pot[!grid->ci][i] = 0.0;
        }
    }

    int steps = compute_potential(grid, threshold);
    printf("Converged after %d\n", steps);    
    compute_density(grid);

    // Print the potential
    char file_name[30] = "test_cond/potential_final.dat";
    FILE* pot_fp = fopen(file_name, "w");
    if (pot_fp == NULL) {
        fprintf(stderr, "Error loading potential output file");
    }
    for (int x = 0; x < grid->width; x++) {
        for (int y = 0; y < grid->height; y++) {
            int index = get_index(grid, x, y);
            fprintf(pot_fp, "%d %d %lf\n", x, y, grid->pot[grid->ci][index]);
        }
    }

    // Print the density
    char den_name[30] = "test_cond/density_final.dat";
    FILE* den_fp = fopen(den_name, "w");
    if (den_fp == NULL) {
        fprintf(stderr, "Error loading density output file");
    }
    for (int x = 0; x < grid->width; x++) {
        for (int y = 0; y < grid->height; y++) {
            int index = get_index(grid, x, y);
            fprintf(den_fp, "%d %d %lf\n", x, y, grid->density[index]);
        }
    }


    // Check if charge sum is close
    double* charges = (double*)calloc(grid->n_conductors, sizeof(double));
    for (int i = 0; i < grid->width * grid->height; i++) {
        if (grid->conductor[i] != 0) {
            charges[grid->conductor[i]-1] += grid->density[i];
        }
    }
    for (int i = 0; i < grid->n_conductors; i++) {
        printf("Original vs computed: (%d) %lf vs %lf\n", i, grid->charges[i], charges[i]);
    }

    free(pots);
    free(charges);
}