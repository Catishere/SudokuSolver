#include <iostream>
#include <ctime>
#include <fstream>
#include <sstream>
#include <cmath>
#include <thread>
#include <vector>
#include <bitset>
#include <cstring>

#define SET_BIT(mask, n) if (n != 0) mask |= 1UL << (n - 1)
#define UNSET_BIT(mask, n) (mask &= ~(1UL << (n)))
#define CHECK_BIT(mask, n) ((mask >> (n)) & 1U)
#define IS_POWER_OF_2(mask) ((mask != 0) && ((mask & (mask - 1)) == 0))
#define IMP(a, b) (~b | a)
#define N_ONES ((1UL << n) - 1)

int n, n_sqrt;

typedef std::vector<std::pair<std::pair<int, int>, std::vector<int>>> cell_info ;

enum CheckType {
    ROW,
    COLUMN,
    MINIGRID
};

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

int** createGrid()
{
    int **grid = new int*[n];
    for (int i = 0; i < n; ++i)
    {
        grid[i] = new int[n]();
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

void availabilityGrid(int **grid, int **mask_grid) {
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
}

uint32_t countSetBits(uint32_t num)
{
    unsigned int count = 0;
    while (num) {
        count += num & 1;
        num >>= 1;
    }
    return count;
}

std::vector<int> getSetBits(int num) {
    std::vector<int> result;
    for (int i = 0; num; ++i) {
        if (~num & 1)
            result.push_back(i + 1);
        num >>= 1;
    }
    return result;
}

int twins(int **mask_grid, int r, int c, CheckType type) {

    int mini_r = r - r % n_sqrt;
    int mini_c = c - c % n_sqrt;
    uint32_t counter = 0;

    switch (type) {
    case ROW:
        for (int i = 0; i < n; ++i) {
            if (mask_grid[r][c] == mask_grid[r][i])
                counter++;
        }
        break;
    case COLUMN:
        for (int i = 0; i < n; ++i) {
            if (mask_grid[r][c] == mask_grid[i][c])
                counter++;
        }
        break;
    case MINIGRID:
        for (int i = mini_r; i < mini_r + n_sqrt; ++i) {
            for (int j = mini_c; j < mini_c + n_sqrt; ++j) {
                if (mask_grid[r][c] == mask_grid[i][j])
                    counter++;
            }
        }
        break;
    }

    if (countSetBits(mask_grid[r][c]) == counter && (counter == 2 || counter == 3))
        return mask_grid[r][c];
    return 0;
}

void loneRanger(int **grid) {
    int **mask_grid = createGrid();
    availabilityGrid(grid, mask_grid);

    for (int r = 0; r < n; ++r) {
        for (int c = 0; c < n; ++c) {
            if (grid[r][c] != 0)
                continue;

            int mask = N_ONES;
            int filter = twins(mask_grid, r, c, ROW);
            int mini_r = r - r % n_sqrt;
            int mini_c = c - c % n_sqrt;

            for (int i = 0; i < n; ++i) {
                toggle_block(mask, mask_grid[r][i], filter);
            }
            mask = ~mask & N_ONES;
            if (IS_POWER_OF_2(mask) && (mask & mask_grid[r][c]) == 0) {
                grid[r][c] = log2i(mask) + 1;
                continue;
            }

            mask = N_ONES;
            filter = twins(mask_grid, r, c, COLUMN);

            for (int i = 0; i < n; ++i) {
                toggle_block(mask, mask_grid[i][c], filter);
            }

            mask = ~mask & N_ONES;
            if (IS_POWER_OF_2(mask) && (mask & mask_grid[r][c]) == 0) {
                grid[r][c] = log2i(mask) + 1;
                continue;
            }

            mask = N_ONES;
            filter = twins(mask_grid, r, c, MINIGRID);

            for (int i = mini_r; i < mini_r + n_sqrt; ++i) {
                for (int j = mini_c; j < mini_c + n_sqrt; ++j) {
                    toggle_block(mask, mask_grid[i][j], filter);
                }
            }

            mask = ~mask & N_ONES;
            if (IS_POWER_OF_2(mask) && (mask & mask_grid[r][c]) == 0) {
                grid[r][c] = log2i(mask) + 1;
                continue;
            }
        }
    }
    destroyGrid(mask_grid);
}

bool isSolved(int **grid) {

    for (int r = 0; r < n; ++r) {
        for (int c = 0; c < n; ++c) {

            uint32_t mask = 0;
            int mini_r = r - r % n_sqrt;
            int mini_c = c - c % n_sqrt;

            for (int i = 0; i < n; ++i) {
                SET_BIT(mask, grid[i][c]);
            }

            if (mask != N_ONES)
                return false;

            mask = 0;

            for (int i = 0; i < n; ++i) {
                SET_BIT(mask, grid[r][i]);
            }

            if (mask != N_ONES)
                return false;


            mask = 0;

            for (int i = mini_r; i < mini_r + n_sqrt; ++i) {
                for (int j = mini_c; j < mini_c + n_sqrt; ++j) {
                    SET_BIT(mask, grid[i][j]);
                }
            }

            if (mask != N_ONES)
                return false;
        }
    }
    return true;
}

void populateGrid(int **grid)
{
    std::stringstream ss;
    ss << "sudoku" << n << ".txt";
    std::ifstream scin(ss.str());
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++) {
            scin >> std::hex >> grid[i][j];
        }
    }
}

void copyGrid(int **dest, int **src) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            dest[i][j] = src[i][j];
        }
    }
}

bool compareGrid(int **a, int **b) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (a[i][j] != b[i][j])
                return false;
        }
    }
    return true;
}

bool bruteForceSearch(int **grid, cell_info &yolo, size_t offset) {

    if (offset >= yolo.size()) {
        if (isSolved(grid)) {
            return true;
        }
        return false;
    }

    for (auto option : yolo[offset].second) {
        grid[yolo[offset].first.first][yolo[offset].first.second] = option;
        bool child_ret = bruteForceSearch(grid, yolo, offset + 1);
        if (child_ret)
            return true;
    }
    return false;
}

void bruteForce(int **grid)
{
    int **mask_grid = createGrid();
    availabilityGrid(grid, mask_grid);

    cell_info yolo;

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++) {
            if (mask_grid[i][j] != 0)
                yolo.push_back({{i, j}, getSetBits(mask_grid[i][j])});
        }
    }

    if (bruteForceSearch(grid, yolo, 0))
        std::cout << "solved" << std::endl;
    else
        std::cout << "fok" << std::endl;
    destroyGrid(mask_grid);
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
    int **grid_new = createGrid();
    populateGrid(grid);

    clock_t start = clock();
    elimination(grid);
    while (!compareGrid(grid, grid_new)) {
        copyGrid(grid_new, grid);
        loneRanger(grid);
        while (!compareGrid(grid, grid_new)) {
            copyGrid(grid_new, grid);
            loneRanger(grid);
        }
        elimination(grid);
    }

    bruteForce(grid);
    clock_t end  = clock();

    print(grid);
    destroyGrid(grid);
    destroyGrid(grid_new);

    float elapsed = (float) (end - start);
    std::cout << elapsed << " ms" << std::endl;
    return 0;
}
