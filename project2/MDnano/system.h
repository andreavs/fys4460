#ifndef SYSTEM_H
#define SYSTEM_H

#include <vector>
#include <armadillo>

#include <string>
#include <math.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>

#include <configreader.h>

class Cell;
class Atom;
class ConfigReader;
class StatisticsCalculator;

class System
{
    friend class Simulation;
    friend class StatisticsCalculator;
public:
    System(ConfigReader *cfgReader);
    int totalCells;
    std::vector<Cell*> cellList;
    double cellSize;
    int totalAtoms;

private:
    int nx;
    int ny;
    int nz;
    int atomsPerGridPoint;

    double volume;
    std::vector<Atom*> atomList;
    void setVelNormal();
    void setVelUniform();
    void setPosFCC();
    void setupCells();

    double tempInKelvin;
    double b;
    double time, temperature, stddev;



    int cellsInXDir;
    int cellsInYDir;
    int cellsInZDir;
    bool createFCC;
    bool readInitialFromFile;
    bool setVelNormalbool;
    bool setVelUniformbool;



    void vmdPrintSystem(std::string filename);
    void placeAtomsInCells();

    std::string experiment; //experiment to start from

    double fastPressure;
    double fastPotentialEnergy;

};

#endif // SYSTEM_H
