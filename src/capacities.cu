#include "lib/solver.cu"

int main(int argc, char *argv[]) {
    if (argc != 4) {
        printf("Usage: %s conv zoom input\n", argv[0]);
        return -1;
    }

    int width, height, channels;
    unsigned char *img = stbi_load(argv[3], &width, &height, &channels, 0);
    if(img == NULL) {
        fprintf(stderr, "Error in loading the image\n");
        exit(1);
    }

    Grid *grid;
    cudaMallocManaged(&grid, sizeof(Grid));
    from_image(grid, img, width, height);
    stbi_image_free(img);

    int zoom = atoi(argv[2]);
    if (zoom > 1) {
        Grid *zoomed;
        cudaMallocManaged(&zoomed, sizeof(Grid));
        copy_zoomed(zoomed, grid, zoom);
        cudaFree(grid);
        
        grid = zoomed;
    }

    compute_capacities(grid, atof(argv[1]));

    // for (int x = 0; x < grid->width; x++) {
    //     for (int y = 0; y < grid->height; y++) {
    //         int index = get_index(grid, x, y);
    //         fprintf(stdout, "%d %d %d\n", x, y, grid->conductor[index]);
    //     }
    // }

    return 0;
}