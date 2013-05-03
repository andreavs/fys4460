#ifndef STATISTICSCALCULATOR_H
#define STATISTICSCALCULATOR_H

#include <iomanip>
#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>

#include "configreader.h"

class System;

class StatisticsCalculator
{
public:
    StatisticsCalculator(System *system, ConfigReader *cfgReader, std::string experiment);
    void printSystemProperties();
    std::ofstream systemFile;
    double energy, kineticEnergy, potentialEnergy, temperature, pressure, displacement, volume;
    void calculateEnergy();
    void calculateTemperature();
    void calculatePressure();
    void calculatePressureFast();
    void calculateDisplacement();
    void calculateVolume();
    void sampleStats();
    System *mySystem;
    bool calculatePressureQuick;
    int numberOfLoopsX;
};

#endif // STATISTICSCALCULATOR_H
