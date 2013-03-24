#ifndef CELL_H
#define CELL_H

#include <armadillo>
#include <vector>
#include <math.h>
#include <list>
#include <deque>
#include <map>

class Atom;


class Cell
{
public:
    Cell(arma::vec3 cellP, double cellS, int cn);
    std::deque<Atom*,int> atomsInCell;
    std::vector<int> neighbourList;
    std::vector<arma::vec3> neighbourVectors;
    int cellIndex;
    double colorIndex;
    arma::vec3 cellPos;
    double cellSize;

};

#endif // CELL_H
