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
    void calculateForces();
    void calculateForcesNull();
    void calculateForcesLJMIC();
    void calculateForcesCellsLJMIC();
    void calculateForcesGravity();
    void calculateForcesThreeBody();
    void calculateForcesFourBody();
    void calculateForcesSixBody();
    void periodicBoundaries();

    bool useLJ;
    bool useGrav;
    bool useThreeBody;
    bool useFourBody;
    bool useSixBody;
    double gravityXComponent;
    double gravityYComponent;
    double gravityZComponent;

    void singlePairForces(double* singlePair, double *singlePairPotential, double *singleParticlePressure);
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

    int cpus;

    int timeStep;

    System *mySystem;
    ConfigReader *cfgReader;
    StatisticsCalculator *myStats;

    double fastPressure;
    double fastPotentialEnergy;


    arma::mat forces;
};

#endif // SIMULATION_H
