#ifndef SYSTEM_H
#define SYSTEM_H

#include "atom.h"
#include "cell.h"
#include <vector>
#include <armadillo>
#include <iomanip>
#include <sstream>
#include <string>
#include <math.h>

class System
{
public:
    System();
    void vmdPrintSystem(std::string filename);
    void runSimulation();

private:
    int nx;
    int ny;
    int nz;
    int atomsPerGridPoint;
    int totalAtoms;
    std::vector<Atom*> atomList;
    void setVelNormal();
    void setPosFCC();
    std::vector<Cell*> cellList;
    double b;
    double mass, dt, F0, E0, T0, sigma, kb, time;
    arma::mat forces;
    void timeEvolve();
    void calculateForcesNull();
    void calculateForcesLJMIC();
    void calculateForcesCellsLJMIC();
    void periodicBoundaries();
    void placeAtomsInCells();
    int noOfTimeSteps;
    int totalCells;
    int cellsInXDir;
    int cellsInYDir;
    int cellsInZDir;
    int cellSize;
};

#endif // SYSTEM_H
