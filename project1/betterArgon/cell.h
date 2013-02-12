#ifndef CELL_H
#define CELL_H

#include <armadillo>
#include <vector>
#include <math.h>
#include "atom.h"
#include <list>


class Cell
{
public:
    Cell(arma::vec3 cellP, double cellS, int cn);
    std::list<Atom*> atomsInCell;
    std::vector<int> neighbourList;
    int cellIndex;

private:
    arma::vec3 cellPos;
    double cellSize;

};

#endif // CELL_H
