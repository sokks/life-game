#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <unistd.h>

class LifeScene {
    int M, N;
    int *cells;

public:
    LifeScene(int M = 10, int N = 10);
    LifeScene(int M, int N, int *cells);
    LifeScene(const LifeScene& ls);
    LifeScene& operator=(const LifeScene& ls);
    ~LifeScene();

    int GetN() const;
    int *GetCells(int row_start, int row_end) const;
    void SetCells(int row_start, int row_end, int*);
    // int *GetCellsPtr(int row_start) const;
    int *GetRow(int i) const;
    void SetRow(int i, const int *row);
    int *GetRowPtr(int i) const;


    int IsAllDead() const;
    void WriteTo(std::string filename) const;
    int Round(int i_start=-1, int i_end=-1);

private:
    int is_alive(int i, int j) const;
    int count_neighs_alive(int i, int j) const;
    int in_row(int i1, int i2) const;
    int in_matrix(int i, int j) const;
};