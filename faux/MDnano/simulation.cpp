#include <armadillo>
#include <list>
#include <omp.h>

#include "simulation.h"
#include "configreader.h"
#include "system.h"
#include "atom.h"
#include "cell.h"
#include "statisticscalculator.h"

#define totalAtoms mySystem->totalAtoms
#define atomList mySystem->atomList
#define nx mySystem->nx
#define ny mySystem->ny
#define nz mySystem->nz
#define b mySystem->b
#define cellsInXDir mySystem->cellsInXDir
#define cellsInYDir mySystem->cellsInYDir
#define cellsInZDir mySystem->cellsInZDir
#define cellList mySystem->cellList



using namespace std;
using namespace arma;

Simulation::Simulation(string filename)
{
    fn = filename;
    string fncfg = "experiments/" + fn + "/" + fn + ".cfg";
    cfgReader = new ConfigReader(fncfg);
    mySystem = new System(cfgReader);
    myStats = new StatisticsCalculator(mySystem, cfgReader);
    absdfadsa = 2;
    try
    {
      noOfTimeSteps = cfgReader->cfg.lookup("noOfTimeSteps");
      calculateStatistics = cfgReader->cfg.lookup("calculateStatistics");
      dt = cfgReader->cfg.lookup("dt");
      useBerendsenThermostat = cfgReader->cfg.lookup("useBerendsenThermostat");
      useAndersenThermostat = cfgReader->cfg.lookup("useAndersenThermostat");
      numberOfThermostatSteps = cfgReader->cfg.lookup("numberOfThermostatSteps");
      thermostatTemperatureKelvin = cfgReader->cfg.lookup("thermostatTemperatureKelvin");
      createMovie = cfgReader->cfg.lookup("createMovie");
      sampleFrequency = cfgReader->cfg.lookup("sampleFrequency");
    }
    catch(const SettingNotFoundException &nfex)
    {
      cerr << "No 'name' setting in configuration file when read from Simulation." << endl;
    }
    forces = zeros<mat>(3,totalAtoms);
    thermostatTemperature = thermostatTemperatureKelvin/119.74;
    cout << "calculating first time forces... ";
    calculateForcesCellsLJMIC();
    cout << "done!" << endl;


}


void Simulation::timeEvolve(){

    for(int i=0; i<totalAtoms;i++){
        atomList[i]->setVel(atomList[i]->getVel() + forces.col(i)*dt/(2*atomList[i]->getMass()));
        atomList[i]->setPos(atomList[i]->getPos() + atomList[i]->getVel()*dt);
        atomList[i]->setRealPos(atomList[i]->getRealPos() + atomList[i]->getVel()*dt);
    }
    mySystem->time += dt;
    periodicBoundaries();
    mySystem->placeAtomsInCells();
    calculateForcesCellsLJMIC();
    for(int i=0; i<totalAtoms;i++){
        atomList[i]->setVel(atomList[i]->getVel() + forces.col(i)*dt/(2*atomList[i]->getMass()));
    }
}


void Simulation::calculateForcesNull(){
    vec3 tempvec = zeros(3);
    for(int i=0; i<totalAtoms;i++){
        forces.col(i) = tempvec;

    }
}

void Simulation::calculateForcesLJMIC(){
    // leonard-jones force using minimal image convention
    calculateForcesNull(); //zero out the forces vector
    vec3 singlePair;
    double singlePairPotential;
    double pressureThread;
    for(int i=0; i<totalAtoms;i++){
        for(int j=i+1;j<totalAtoms;j++){
            singlePair = (atomList[i]->getPos() - atomList[j]->getPos());
            singlePairForces(singlePair, &singlePairPotential, &pressureThread);
            forces.col(i) += singlePair;
            forces.col(j) -= singlePair;
        }
    }
}

void Simulation::singlePairForces(vec3& singlePair, double *singlePairPotential, double *pressureThread){
    double r2inv;
    double r6inv;
    double r8inv;
    double factor;
//    for (int j=0; j<3; j++){
//        singlePair(j) = (abs(singlePair(j)) < abs(singlePair(j) + nx*b))*(abs(singlePair(j)) < abs(singlePair(j) - nx*b))*singlePair(j)
//                + (abs(singlePair(j) + nx*b) < abs(singlePair(j)))*(abs(singlePair(j) + nx*b) < abs(singlePair(j) - nx*b))*(singlePair(j)+nx*b)
//                + (abs(singlePair(j) - nx*b) < abs(singlePair(j)))*(abs(singlePair(j) - nx*b) < abs(singlePair(j) + nx*b))*(singlePair(j)-nx*b);
//    }

    r2inv = 1./(singlePair(0)*singlePair(0) +  singlePair(1)*singlePair(1) + singlePair(2)*singlePair(2));
    //r2 = singlePair(0)*singlePair(0) +  singlePair(1)*singlePair(1) + singlePair(2)*singlePair(2);
    //r2 = max(dot(singlePair,singlePair),0.81);
    r6inv = r2inv*r2inv*r2inv;
    r8inv = r2inv*r6inv;
    factor = (24.*r8inv)*(2*r6inv - 1);
    *singlePairPotential = 4.0*r6inv*(r6inv - 1.0);
    *pressureThread += factor*( singlePair(0)*singlePair(0) + singlePair(1)*singlePair(1) + singlePair(2)*singlePair(2) );
    singlePair(0) = factor*singlePair(0);
    singlePair(1) = factor*singlePair(1);
    singlePair(2) = factor*singlePair(2);
    //singlePair = singlePair*(24./r8)*(2/r6 - 1);
}

void Simulation::runSimulation(){
    ostringstream number;
    // "experiments/" + fn + "/" + fn + ".cfg";
    string filename = "experiments/" + fn + "/results/results";
    string saveFilename;
    string extension = ".xyz";
    //mySystem->vmdPrintSystem(saveFilename);
    for(int i=0; i<noOfTimeSteps;i++){
        calculateForcesNull();
        cout << "solving for time step: " << i << " out of: " << noOfTimeSteps << "...";
        //some jargon to make file extensions 000, 001, 002, etc.
        filename = "experiments/" + fn + "/results/results";
        number << setw(4) << setfill('0') << i;
        saveFilename = number.str();
        number.str(""); //TODO: This is bad use of ostringstream, not meant to be reused!
        saveFilename = filename.append(saveFilename).append(extension);
        //jargon done!
        if(createMovie){
            mySystem->vmdPrintSystem(saveFilename);
        }
        if(i==noOfTimeSteps-1){
            saveFilename = "lastState.xyz";
            mySystem->vmdPrintSystem(saveFilename);
        }
        if(calculateStatistics && (i%sampleFrequency == 0)){
            myStats->sampleStats();
            myStats->printSystemProperties();
        }
        fastPressure = 0.0;
        fastPotentialEnergy = 0.0;
        //calculateForcesCellsLJMIC();
        timeEvolve();
        if(useAndersenThermostat && i < numberOfThermostatSteps){
            andersenThermostat(thermostatTemperature);
        }
        else if(useBerendsenThermostat && i < numberOfThermostatSteps){
            berendsenThermostat(thermostatTemperature);
        }
        cout << " done!" << endl;
    }
    //some jargon to make file extensions 000, 001, 002, etc.
    filename = "results";
    number << setw(4) << setfill('0') << noOfTimeSteps;
    saveFilename = number.str();
    number.str(""); //TODO: This is bad use of ostringstream, not meant to be reused!
    saveFilename = filename.append(saveFilename).append(extension);
    //jargon done!
    mySystem->vmdPrintSystem(saveFilename);
    //calculate kinetic energy

    myStats->sampleStats();
    myStats->printSystemProperties();
    myStats->systemFile.close();
}

void Simulation::periodicBoundaries(){
    vec3 modulusVec;
    double xtest = b*nx;
    double ytest = b*ny;
    double ztest = b*nz;
    for(int i=0; i<totalAtoms;i++){
        modulusVec = atomList[i]->getPos();
        modulusVec(0) = fmod(modulusVec(0),xtest);
        if(modulusVec(0) < 0){modulusVec(0) += xtest;}
        modulusVec(1) = fmod(modulusVec(1),ytest);
        if(modulusVec(1) < 0){modulusVec(1) += ytest;}
        modulusVec(2) = fmod(modulusVec(2),ztest);
        if(modulusVec(2) < 0){modulusVec(2) += ztest;}
        atomList[i]->setPos(modulusVec);
    }
}



void Simulation::calculateForcesCellsLJMIC(){
    // leonard-jones force using minimal image convention
    calculateForcesNull(); //zero out the forces vector
#pragma omp parallel num_threads(1)
    {
    vec3 singlePair;
    int i;
    int thisIndex;
    int otherIndex;
    int neighbourPos;
    mat forcesThread = zeros(3,totalAtoms);
    double potentialEnergyThread = 0.0;
    double singlePairPotential;
    double pressureThread;
//#pragma omp parallel private(singlePair,singlePairPotential,pressureThread,i,thisIndex,otherIndex,dummyk,jint,j,k,neighbourPos,counter,forcesThread,potentialEnergyThread) num_threads(4)
//    {
//        forcesThread = zeros(3,totalAtoms);
//        potentialEnergyThread = 0.0;
//        pressureThread = 0.0;
    #pragma omp for
        for(int zCellPos=0; zCellPos<cellsInXDir;zCellPos++){
            for(int yCellPos=0; yCellPos<cellsInYDir;yCellPos++){
                for(int xCellPos=0; xCellPos<cellsInYDir;xCellPos++){
                    i = zCellPos*cellsInXDir*cellsInYDir + yCellPos*cellsInXDir + xCellPos;
                    for(auto j=cellList[i]->atomsInCell.begin(); j!= cellList[i]->atomsInCell.end(); ++j){

                        //particles in own cell
                        thisIndex = (*j)->getSystemIndex();
                        auto dummyk = j;
                        advance(dummyk,1);
                        for(auto k=dummyk; k != cellList[i]->atomsInCell.end(); ++k){
                            otherIndex = (*k)->getSystemIndex();
                            singlePair(0) = (atomList[thisIndex]->getPos()(0) - atomList[otherIndex]->getPos()(0));
                            singlePair(1) = (atomList[thisIndex]->getPos()(1) - atomList[otherIndex]->getPos()(1));
                            singlePair(2) = (atomList[thisIndex]->getPos()(2) - atomList[otherIndex]->getPos()(2));

                            singlePairForces(singlePair,&singlePairPotential, &pressureThread);
                            potentialEnergyThread += singlePairPotential;
                            forcesThread(0,thisIndex) += singlePair(0); forcesThread(0,otherIndex) -= singlePair(0);
                            forcesThread(1,thisIndex) += singlePair(1); forcesThread(1,otherIndex) -= singlePair(1);
                            forcesThread(2,thisIndex) += singlePair(2); forcesThread(2,otherIndex) -= singlePair(2);
                        }

                        //particles in neighbouring cells
                        //each cell points to 13 neighbouring cells, to avoid double counting
                        //the cell points to cells which have (x,y,z)-indices increased by 1 (mod N)
                        //each cell is initialized with a neighbour vector<int> called neighbourList

                        for(int neighbour = 0; neighbour<cellList[i]->neighbourList.size(); neighbour++){
                            neighbourPos = cellList[i]->neighbourList[neighbour];
                            for(auto k=cellList[neighbourPos]->atomsInCell.begin(); k!= cellList[neighbourPos]->atomsInCell.end(); ++k){
                                otherIndex = (*k)->getSystemIndex();
                                singlePair(0) = (atomList[thisIndex]->getPos()(0) - atomList[otherIndex]->getPos()(0)) - cellList[i]->neighbourVectors[neighbour](0);
                                singlePair(1) = (atomList[thisIndex]->getPos()(1) - atomList[otherIndex]->getPos()(1)) - cellList[i]->neighbourVectors[neighbour](1);
                                singlePair(2) = (atomList[thisIndex]->getPos()(2) - atomList[otherIndex]->getPos()(2)) - cellList[i]->neighbourVectors[neighbour](2);

                                //singlePair = (atomList[thisIndex]->getPos() - atomList[otherIndex]->getPos()) - cellList[i]->neighbourVectors[neighbour];
                                singlePairForces(singlePair, &singlePairPotential, &pressureThread);
                                potentialEnergyThread += singlePairPotential;
                                forcesThread(0,thisIndex) += singlePair(0); forcesThread(0,otherIndex) -= singlePair(0);
                                forcesThread(1,thisIndex) += singlePair(1); forcesThread(1,otherIndex) -= singlePair(1);
                                forcesThread(2,thisIndex) += singlePair(2); forcesThread(2,otherIndex) -= singlePair(2);
                            }
                        }
                    }
                }
            }
        }
        #pragma omp critical
        {
            fastPotentialEnergy = fastPotentialEnergy + potentialEnergyThread;
            fastPressure = fastPressure + pressureThread;
            for(int i=0;i<3;i++){
                for(int j=0; j<totalAtoms;j++){
                    forces(i,j) = forces(i,j) + forcesThread(i,j);
                }
            }
        }
    }
    myStats->pressure = fastPressure;
    myStats->potentialEnergy = fastPotentialEnergy;
}


void Simulation::berendsenThermostat(double temp){
    double tau = 10*dt;
    double gamma = sqrt(1 + dt/tau*(temp/mySystem->temperature - 1));
    for(int i=0; i<totalAtoms; i++){
        atomList[i]->setVel(gamma*atomList[i]->getVel());
    }
}

void Simulation::andersenThermostat(double temp){
    vec3 newVel = zeros(3);
    vec randnumber = randu(totalAtoms);
    double tau = 10*dt;
    double test = dt/tau;
    for(int i=0; i<totalAtoms;i++){
        if(randnumber[i] < test){
            newVel = sqrt(temp)*randn(3);
            atomList[i]->setVel(newVel);
        }

    }
}


