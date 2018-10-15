#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <unistd.h>

#include "life_scene.h"

LifeScene::LifeScene(int M, int N)
{
    this->M = M;
    this->N = N;
    cells = new int[M * N];
    for (int i = 0; i < M * N; i++)
    {
        cells[i] = std::rand() % 2;
    }
}

LifeScene::LifeScene(int M, int N, int *cells) {
    this->M = M;
    this->N = N;
    this->cells = new int[M * N];
    for (int i = 0; i < M * N; i++)
    {
        this->cells[i] = cells[i];
    }
}

LifeScene::LifeScene(const LifeScene &ls)
{
    M = ls.M;
    N = ls.N;
    cells = new int[M * N];
    for (int i = 0; i < M * N; i++)
    {
        cells[i] = ls.cells[i];
    }
}

LifeScene &LifeScene::operator=(const LifeScene &ls)
{
    M = ls.M;
    N = ls.N;
    delete[] cells;
    cells = new int[M * N];
    for (int i = 0; i < M * N; i++)
    {
        cells[i] = ls.cells[i];
    }
    return *this;
}

LifeScene::~LifeScene()
{
    delete[] cells;
}

int LifeScene::IsAllDead() const
{
    for (int i = 0; i < M * N; i++)
    {
        if (cells[i])
        {
            return 0;
        }
    }
    return 1;
}

void LifeScene::WriteTo(std::string filename) const
{
    std::ofstream fout(filename, std::ios::trunc | std::ios::out);
    // fout << N << " " << M << std::endl;
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < M; j++)
        {
            fout << cells[i * M + j] << " ";
        }
        fout << std::endl;
    }
    fout.close();
}

int LifeScene::Round(int i_start, int i_end)
{
    if (i_start == -1) {
        i_start = 0;
    }
    if (i_end == -1) {
        i_end = N;
    }

    int *new_cells = new int[M * N];
    for (int i = i_start; i < i_end; i++)
    {
        for (int j = 0; j < M; j++)
        {
            int neighs_alive = this->count_neighs_alive(i, j);
            if (is_alive(i, j))
            {
                new_cells[i * M + j] = ((neighs_alive == 2) || (neighs_alive == 3)) ? 1 : 0;
            }
            else
            {
                new_cells[i * M + j] = (neighs_alive == 3) ? 1 : 0;
            }
        }
    }
    int all_the_same = 1;
    for (int i = i_start*M; i < M * (i_end); i++)
    {
        if (cells[i] != new_cells[i])
        {
            all_the_same = 0;
            break;
        }
    }
    if (IsAllDead() || all_the_same)
    {
        return 1;
    }
    delete[] cells;
    cells = new_cells;
    return 0;
}

int LifeScene::is_alive(int i, int j) const
{
    return cells[i * M + j];
}

int LifeScene::count_neighs_alive(int i, int j) const
{
    int i_left = i - 1;
    int i_right = i + 1;
    int i_down = i + M;
    int i_up = i - M;

    int i_upleft = i - M - 1;
    int i_upright = i - M + 1;
    int i_downleft = i + M - 1;
    int i_downright = i + M + 1;

    int n_alive = 0;
    int range[] = {-1, 0, 1};
    for (int k = 0; k < 3; k++)
    {
        for (int l = 0; l < 3; l++)
        {
            int delta_i = range[k];
            int delta_j = range[l];
            if (!delta_i && !delta_j)
            {
                continue;
            }
            n_alive += in_matrix(i + delta_i, j + delta_j) && is_alive(i + delta_i, j + delta_j);
        }
    }

    return n_alive;
}

int LifeScene::in_row(int i1, int i2) const
{
    return (i1 / M == i2 / M) ? 1 : 0;
}

int LifeScene::in_matrix(int i, int j) const
{
    return ((i >= 0) && (i < N) && (j >= 0) && (j < M)) ? 1 : 0;
}

int *LifeScene::GetCells(int row_start, int row_end) const {
    int *block = new int[M * (row_end - row_start)];
    int k = 0;
    for (int i = row_start; i < row_end; i++) {
        for (int j = 0; j < M; j++) {
            block[k] = cells[i*M + j];
            k++;
        }
    }

    return block;
}

int *LifeScene::GetRow(int i) const {
    int *row = new int[M];
    for (int j = 0; j < M; j++) {
        row[j] = cells[i*M + j];
    }

    return row;
}

void LifeScene::SetRow(int i, const int *row) {
    for (int j = 0; j < M; j++) {
        cells[i*M + j] = row[j];
    }
}

void LifeScene::SetCells(int row_start, int row_end, int *rows) {
    int k = 0;
    for (int i = row_start; i < row_end; i++) {
        for (int j = 0; j < M; j++) {
            cells[i*M+j] = rows[k];
            k++;
        }
    }
}

int *LifeScene::GetRowPtr(int i) const {
    return cells+i*M;
}
