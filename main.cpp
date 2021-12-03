#include <iostream>
#include <ctime>
#include <fstream>
#include <sstream>
#include <cmath>
#include <thread>
#include <mutex>
#include <vector>
#include <bitset>
#include <cstring>
#include <iomanip>

#define SET_BIT(mask, n) if (n != 0) mask |= 1UL << (n - 1)
#define UNSET_BIT(mask, n) (mask &= ~(1UL << (n)))
#define CHECK_BIT(mask, n) ((mask >> (n)) & 1U)
#define IS_POWER_OF_2(mask) ((mask != 0) && ((mask & (mask - 1)) == 0))
#define IMP(a, b) (~b | a)
#define N_ONES ((1UL << n) - 1)

int n, n_sqrt;

enum CheckType {
    ROW,
    COLUMN,
    MINIGRID
};

void print(uint32_t **arr)
{
    for (int i = 0; i < n; i++)
    {
        if (i % n_sqrt == 0)
            std::cout << std::string(n * 2 + 6, '-') << std::endl;
        for (int j = 0; j < n; j++) {
            if (j % n_sqrt == 0)
                std::cout << "| ";
            std::cout << std::setw(2) <<  std::left << arr[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

static inline int log2i(int x) {
    return sizeof(int) * 8 - __builtin_clz(x) - 1;
}

uint32_t** createGrid()
{
    uint32_t **grid = new uint32_t*[n];
    for (int i = 0; i < n; ++i)
    {
        grid[i] = new uint32_t[n]();
    }
    return grid;
}

void destroyGrid(uint32_t **grid)
{
    for (int i = 0; i < n; ++i)
    {
        delete [] grid[i];
    }
    delete [] grid;
}

void elimination(uint32_t **grid) {
    for (int r = 0; r < n; ++r) {
        for (int c = 0; c < n; ++c) {
            if (grid[r][c] != 0)
                continue;

            uint32_t mask = N_ONES;
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

static inline void toggle_block(uint32_t &mask, uint32_t grid_mask, uint32_t &filter) {
    if (grid_mask == 0)
        return;
    uint32_t new_mask = ~(~mask ^ ~grid_mask); // xnor stack new, exclude repetitives
    filter |= ~mask & ~grid_mask; // 1 if bit is excluded
    filter &= N_ONES;
    mask = new_mask | filter;
    mask &= N_ONES;
}

void availabilityGrid(uint32_t **grid, uint32_t **mask_grid) {
    for (int r = 0; r < n; ++r) {
        for (int c = 0; c < n; ++c) {
            if (grid[r][c] != 0)
                continue;

            uint32_t mask = 0;
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

int countSetBits(int num)
{
    int count = 0;
    while (num) {
        count += num & 1;
        num >>= 1;
    }
    return count;
}

std::vector<int> getSetBits(uint32_t num) {
    std::vector<int> result;
    for (int i = 0; num; ++i) {
        if (~num & 1)
            result.push_back(i + 1);
        num >>= 1;
    }
    return result;
}

uint32_t twins(uint32_t **mask_grid, int r, int c, CheckType type) {

    int mini_r = r - r % n_sqrt;
    int mini_c = c - c % n_sqrt;
    int counter = 0;

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

    if (countSetBits(mask_grid[r][c]) == counter
            && counter > 1
            && counter < n)
        return mask_grid[r][c];
    return 0;
}

void loneRanger(uint32_t **grid) {
    uint32_t **mask_grid = createGrid();
    availabilityGrid(grid, mask_grid);

    for (int r = 0; r < n; ++r) {
        for (int c = 0; c < n; ++c) {
            if (grid[r][c] != 0)
                continue;

            uint32_t mask = N_ONES;
            uint32_t filter = twins(mask_grid, r, c, ROW);
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

void populateGrid(uint32_t **grid)
{
    std::stringstream ss;
    ss << "sudoku" << n << ".txt";
    std::ifstream scin(ss.str());
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++) {
            scin >> grid[i][j];
        }
    }
}

void copyGrid(uint32_t **dest, uint32_t **src) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            dest[i][j] = src[i][j];
        }
    }
}

bool compareGrid(uint32_t **a, uint32_t **b) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (a[i][j] != b[i][j])
                return false;
        }
    }
    return true;
}

bool isSafe(uint32_t **grid, int row, int col, uint32_t num)
{
    for (int x = 0; x <= n - 1; x++)
        if (grid[row][x] == num)
            return false;

    for (int x = 0; x <= n - 1; x++)
        if (grid[x][col] == num)
            return false;

    int startRow = row - row % n_sqrt;
    int startCol = col - col % n_sqrt;

    for (int i = 0; i < n_sqrt; i++)
        for (int j = 0; j < n_sqrt; j++)
            if (grid[i + startRow][j + startCol] == num)
                return false;

    return true;
}

bool solveSuduko(uint32_t **grid, int row, int col)
{
    if (row == n - 1 && col == n)
        return true;

    if (col == n) {
        row++;
        col = 0;
    }
    if (grid[row][col] > 0)
        return solveSuduko(grid, row, col + 1);

    for (int num = 1; num <= n; num++)
    {
        if (isSafe(grid, row, col, num))
        {
            grid[row][col] = num;

            if (solveSuduko(grid, row, col + 1))
                return true;
        }
        grid[row][col] = 0;
    }
    return false;
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

    uint32_t **grid = createGrid();
    uint32_t **grid_new = createGrid();
    populateGrid(grid);

    clock_t start = clock();
    elimination(grid);
    while (!compareGrid(grid, grid_new)) {
        copyGrid(grid_new, grid);
        loneRanger(grid);
    }

    auto ret = solveSuduko(grid, 0, 0);

    clock_t end  = clock();

    if (!ret)
        std::cout << "no solution  exists " << std::endl;

    print(grid);
    destroyGrid(grid);
    destroyGrid(grid_new);

    float elapsed = (float) (end - start);
    std::cout << elapsed << " ms" << std::endl;
    return 0;
}
