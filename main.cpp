#include <iostream>
#include <ctime>
#include <fstream>
#include <sstream>
#include <cmath>
#include <thread>
#include <bitset>
#include <cstring>

#define SET_BIT(mask, n) if (n != 0) mask |= 1UL << (n - 1)
#define UNSET_BIT(mask, n) (mask &= ~(1UL << (n)))
#define CHECK_BIT(mask, n) ((mask >> (n)) & 1U)
#define IS_POWER_OF_2(mask) ((mask != 0) && ((mask & (mask - 1)) == 0))
#define IMP(a, b) (~b | a)
#define N_ONES ((1UL << n) - 1)

int n, n_sqrt;

void print(int **arr)
{
    for (int i = 0; i < n; i++)
    {
        if (i % n_sqrt == 0)
            std::cout << std::string(n * 2 + 6, '-') << std::endl;
        for (int j = 0; j < n; j++) {
            if (j % n_sqrt == 0)
                std::cout << "| ";
            std::cout << arr[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

static inline int log2i(int x) {
    return sizeof(int) * 8 - __builtin_clz(x) - 1;
}

void elimination(int **grid) {
    for (int r = 0; r < n; ++r) {
        for (int c = 0; c < n; ++c) {
            if (grid[r][c] != 0)
                continue;

            int mask = N_ONES;
            int mini_r = r - r % n_sqrt;
            int mini_c = c - c % n_sqrt;

            for (int i = 0; i < n; ++i) {
                UNSET_BIT(mask, grid[i][c]);
            }

            for (int i = 0; i < n; ++i) {
                UNSET_BIT(mask, grid[r][i]);
            }

            for (int i = mini_r; i < n_sqrt; ++i) {
                for (int j = mini_c; j < n_sqrt; ++j) {
                    UNSET_BIT(mask, grid[i][j]);
                }
            }
            // check if mask is a power of 2 and get the index with log2i
            if (IS_POWER_OF_2(mask))
                grid[r][c] = log2i(mask);
        }
    }
}

static inline void toggle_block(int &mask, int grid_mask, int &filter) {
    if (grid_mask == 0)
        return;
    int new_mask = ~(~mask ^ ~grid_mask); // xnor stack new, exclude repetitives
    filter |= ~mask & ~grid_mask; // 1 if bit is excluded
    filter &= N_ONES;
    mask = new_mask | filter;
    mask &= N_ONES;
}

void loneRanger(int **grid) {
    int mask_grid[n][n];
    memset(mask_grid, 0, sizeof(mask_grid));
    for (int r = 0; r < n; ++r) {
        for (int c = 0; c < n; ++c) {
            if (grid[r][c] != 0)
                continue;

            int mask = 0;
            int mini_r = r - r % n_sqrt;
            int mini_c = c - c % n_sqrt;

            for (int i = 0; i < n; ++i) {
                SET_BIT(mask, grid[i][c]);
            }

            for (int i = 0; i < n; ++i) {
                SET_BIT(mask, grid[r][i]);
            }

            for (int i = mini_r; i < mini_r + n_sqrt; ++i) {
                for (int j = mini_c; j < mini_c + n_sqrt; ++j) {
                    SET_BIT(mask, grid[i][j]);
                }
            }
            mask_grid[r][c] = mask;
        }
    }

    for (int r = 0; r < n; ++r) {
        for (int c = 0; c < n; ++c) {
            if (grid[r][c] != 0)
                continue;

            int mask = N_ONES;
            int filter = 0;
            int mini_r = r - r % n_sqrt;
            int mini_c = c - c % n_sqrt;

            for (int i = 0; i < n; ++i) {
                toggle_block(mask, mask_grid[r][i], filter);
            }
            mask = ~mask & N_ONES;
            if (IS_POWER_OF_2(mask)) {
                grid[r][c] = log2i(mask) + 1;
                continue;
            }

            mask = N_ONES;
            filter = 0;

            for (int i = 0; i < n; ++i) {
                toggle_block(mask, mask_grid[i][c], filter);
            }

            mask = ~mask & N_ONES;
            if (IS_POWER_OF_2(mask)) {
                grid[r][c] = log2i(mask) + 1;
                continue;
            }

            mask = N_ONES;
            filter = 0;

            for (int i = mini_r; i < mini_r + n_sqrt; ++i) {
                for (int j = mini_c; j < mini_c + n_sqrt; ++j) {
                    toggle_block(mask, mask_grid[i][j], filter);
                }
            }

            mask = ~mask & N_ONES;
            if (IS_POWER_OF_2(mask)) {
                grid[r][c] = log2i(mask) + 1;
                continue;
            }
        }
    }
}

int** createGrid()
{
    int **grid = new int*[n];
    for (int i = 0; i < n; ++i)
    {
        grid[i] = new int[n];
    }
    return grid;
}

void destroyGrid(int **grid)
{
    for (int i = 0; i < n; ++i)
    {
        delete [] grid[i];
    }
    delete [] grid;
}

void populateGrid(int **grid)
{
    std::stringstream ss;
    ss << "sudoku" << n << "ranger.txt";
    std::ifstream scin(ss.str());
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++) {
            scin >> std::hex >> grid[i][j];
        }
    }
}

int main(int argc, char **argv)
{
    if (argc < 2)
        return 1;

    n = std::stoi(argv[1]);

    if (n <= 0) return 1;

    n_sqrt = std::sqrt(n);

    if (n_sqrt * n_sqrt != n) {
        std::cout << "N should be a square root" << std::endl;
        return 1;
    }

    int **grid = createGrid();
    populateGrid(grid);

    clock_t start = clock();
    //elimination(grid);
    loneRanger(grid);
    clock_t end  = clock();

    print(grid);
    destroyGrid(grid);

    float elapsed = (float) (end - start);
    std::cout << elapsed << " ms" << std::endl;
    return 0;
}
