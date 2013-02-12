#include "cell.h"
using namespace arma;

Cell::Cell(vec3 cellP, double cellS, int cn)
{
    cellPos = cellP;
    cellSize = cellS;
    cellIndex = cn;

}
