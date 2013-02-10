#ifndef CELL_H
#define CELL_H

#include <armadillo>
#include <vector>
#include <math.h>
#include "atom.h"


class Cell
{
public:
    Cell(arma::vec3 cellP, double cellS);
    std::vector<Atom*> atomsInCell;

private:
    arma::vec3 cellPos;
    double cellSize;

};

#endif // CELL_H
