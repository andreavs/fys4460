#ifndef CELL_H
#define CELL_H

#include <armadillo>
#include "atom.h"
#include <vector>
#include <math.h>


class Cell
{
public:
    Cell(double cellP, double cellS);
    std::vector<Atom*> atomsInCell;

private:
    arma::vec3 cellPos;
    double cellSize;

};

#endif // CELL_H
