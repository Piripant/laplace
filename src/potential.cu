#include "lib/solver.cu"

int main(int argc, char *argv[]) {
    if (argc != 5) {
        printf("Usage: %s threshold i_pot o_pot o_charges\n", argv[0]);
        return -1;
    }

    int width, height, channels;
    unsigned char *img = stbi_load(argv[2], &width, &height, &channels, 0);
    if(img == NULL) {
        fprintf(stderr, "Error in loading the image\n");
        exit(1);
    }

    Grid *grid;
    cudaMallocManaged(&grid, sizeof(Grid));
    init_grid(grid, width, height);

    size_t img_size = width * height * channels;
    for(int index = 0; index < img_size; index += channels) {
        int color = img[index];
        if (color != 0) {
            grid->conductor[index] = 1;
            grid->pot[grid->ci][index] = (double)color / 256.0;
            grid->pot[!grid->ci][index] = (double)color / 256.0;
        }
    }
    stbi_image_free(img);

    // Compute potential
    double threshold = atof(argv[1]);
    int count = compute_potential(grid, threshold);

    printf("Converged after %d steps\n", count);

    // Print the potential
    FILE* pot_fp = fopen(argv[3], "w");
    if (pot_fp == NULL) {
        fprintf(stderr, "Error loading potential output file");
        return -1;
    }
    for (int x = 0; x < grid->width; x++) {
        for (int y = 0; y < grid->height; y++) {
            int index = get_index(grid, x, y);
            fprintf(pot_fp, "%d %d %lf\n", x, y, grid->pot[grid->ci][index]);
        }
    }

    // Compute the distribution of charges
    compute_density(grid);

    printf("Priting charges\n");
    // Print the charges
    FILE* charge_fp = fopen(argv[4], "w");
    if (charge_fp == NULL) {
        fprintf(stderr, "Error loading charges output file");
        return -1;
    }
    
    for (int x = 0; x < grid->width; x++) {
        for (int y = 0; y < grid->height; y++) {
            int index = get_index(grid, x, y);
            fprintf(charge_fp, "%d %d %lf\n", x, y, grid->density[index]);
            if (grid->density[index] != 0.0) {
                printf("Not 0\n");
            }
        }
    }
    
    return 0;
}