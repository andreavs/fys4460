#ifndef SIMULATION_H
#define SIMULATION_H

#include <string>
#include <armadillo>

class System;
class ConfigReader;
class StatisticsCalculator;

class Simulation
{
public:
    Simulation(std::string filename);
    void vmdPrintSystem(std::string filename);
    void runSimulation();

    void timeEvolve();
    void calculateForcesNull();
    void calculateForcesLJMIC();
    void calculateForcesCellsLJMIC();
    void periodicBoundaries();

    void singlePairForces(double* singlePair, double *singlePairPotential, double *pressureThread);
    void berendsenThermostat(double temp);
    void andersenThermostat(double temp);
    double temperatureBath();

    bool calculateStatistics;
    int noOfTimeSteps;
    double dt;
    bool useBerendsenThermostat;
    bool useAndersenThermostat;
    int numberOfThermostatSteps;
    double thermostatTemperatureKelvin;
    double thermostatTemperature;
    std::string xyzFile;
    bool createMovie;

    std::string fn; //config file to be read
    int sampleFrequency; //how often to sample stats (in time step)
    bool sampleNow;

    System *mySystem;
    ConfigReader *cfgReader;
    StatisticsCalculator *myStats;

    double fastPressure;
    double fastPotentialEnergy;

    arma::mat forces;
};

#endif // SIMULATION_H
