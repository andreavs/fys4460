#include "statisticscalculator.h"
#include "system.h"
#include "atom.h"

using namespace std;
using namespace arma;

StatisticsCalculator::StatisticsCalculator(System *system, ConfigReader* cfgReader)
{
    try
    {
        calculatePressureQuick = cfgReader->cfg.lookup("calculatePressureQuick");

    }
    catch(const SettingNotFoundException &nfex)
    {
      cerr << "No 'name' setting in configuration file when read from StatisticsCalculator." << endl;
    }
    string filename = "systemresults.txt";
    systemFile.open(filename.c_str());
    systemFile << "(time) (Energy) (Kinetic Energy) (Potential Energy) (Pressure) (Temperature) " << std::endl;
    pressure = 0.0;
    mySystem = system;
}


void StatisticsCalculator::printSystemProperties(){
    //this function saves the properties of the system such as energy pressure etc
    systemFile << mySystem->time << " " << energy << " " << kineticEnergy << " " << potentialEnergy
               << " " << pressure << " " << temperature << " " << displacement << " " << std::endl;
}


void StatisticsCalculator::calculateEnergy(){
    kineticEnergy = 0;
    for(int i=0; i<mySystem->totalAtoms; i++){
        kineticEnergy += 0.5*mySystem->atomList[i]->getMass()*dot(mySystem->atomList[i]->getVel(), mySystem->atomList[i]->getVel());
    }
    //potentialEnergy = mySystem->fastPotentialEnergy;

    energy = kineticEnergy + potentialEnergy;
}

void StatisticsCalculator::calculateTemperature(){
    temperature = 2./(3.*mySystem->totalAtoms)*kineticEnergy;
}

void StatisticsCalculator::calculatePressureFast(){
    //pressure = mySystem->fastPressure;
    pressure = mySystem->totalAtoms/volume*temperature + (1./(3.*volume))*pressure;
}

void StatisticsCalculator::calculatePressure(){

}

void StatisticsCalculator::calculateDisplacement(){
    displacement = 0;
    for(int i=0; i<mySystem->totalAtoms; i++){
        displacement += dot(mySystem->atomList[i]->getRealPos() - mySystem->atomList[i]->getInitialPos(), mySystem->atomList[i]->getRealPos() - mySystem->atomList[i]->getInitialPos());
    }
    displacement /= mySystem->totalAtoms;
}

void StatisticsCalculator::calculateVolume(){

}

void StatisticsCalculator::sampleStats(){
    calculateEnergy();
    calculateTemperature();
    calculatePressure();
    calculateDisplacement();
    calculateVolume();
}
