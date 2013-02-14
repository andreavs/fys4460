#ifndef SYSTEM_H
#define SYSTEM_H


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
    int totalCells;
    std::vector<Cell*> cellList;
    double cellSize;

private:
    int nx;
    int ny;
    int nz;
    int atomsPerGridPoint;
    int totalAtoms;
    std::vector<Atom*> atomList;
    void setVelNormal();
    void setPosFCC();
    double b;
    double mass, dt, F0, E0, T0, sigma, kb, time;
    arma::mat forces;
    void timeEvolve();
    void calculateForcesNull();
    void calculateForcesLJMIC();
    void calculateForcesCellsLJMIC();
    void periodicBoundaries();
    void placeAtomsInCells();
    void singlePairForces(arma::vec3 &singlePair);
    int noOfTimeSteps;
    int cellsInXDir;
    int cellsInYDir;
    int cellsInZDir;

};

#endif // SYSTEM_H
