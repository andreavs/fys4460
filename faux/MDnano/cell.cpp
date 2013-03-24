#include "cell.h"
using namespace arma;

Cell::Cell(vec3 cellP, double cellS, int cn)
{
    cellPos = cellP;
    cellSize = cellS;
    cellIndex = cn;
    vec v1 = randu(1);
    colorIndex = v1[0];

}
