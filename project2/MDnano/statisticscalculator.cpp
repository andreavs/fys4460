#include "statisticscalculator.h"
#include "system.h"
#include "atom.h"

using namespace std;
using namespace arma;

StatisticsCalculator::StatisticsCalculator(System *system, ConfigReader* cfgReader, string experiment)
{
    string filename;
    try
    {
        calculatePressureQuick = cfgReader->cfg.lookup("calculatePressureQuick");
        //string fn = cfgReader->cfg.lookup("experiment");
        //filename = fn;

    }
    catch(const SettingNotFoundException &nfex)
    {
      cerr << "No 'name' setting in configuration file when read from StatisticsCalculator." << endl;
    }
    filename = "experiments/" + experiment + "/results/systemresults.txt";
    //string filename = "systemresults.txt";
    systemFile.open(filename.c_str());
    systemFile << "(time) (Energy) (Kinetic Energy) (Potential Energy) (Pressure) (Temperature) " << std::endl;
    pressure = 0.0;
    mySystem = system;
    volume = pow(mySystem->b,3)*mySystem->nx*mySystem->ny*mySystem->nz;
    numberOfLoopsX = 0;
}


void StatisticsCalculator::printSystemProperties(){
    //this function saves the properties of the system such as energy pressure etc
    //cout << potentialEnergy << endl;
    systemFile << mySystem->time << " " << energy << " " << kineticEnergy << " " << potentialEnergy
               << " " << pressure << " " << temperature << " " << displacement << " " << numberOfLoopsX << " " << std::endl;
}


void StatisticsCalculator::calculateEnergy(){
    kineticEnergy = 0;

    for(int i=0; i<mySystem->totalAtoms; i++){
        if(!mySystem->atomList[i]->isFrozen){
            kineticEnergy += 0.5*mySystem->atomList[i]->getMass()*dot(mySystem->atomList[i]->getVel(), mySystem->atomList[i]->getVel());
        }
    }
    //potentialEnergy = mySystem->fastPotentialEnergy;
    //cout << potentialEnergy << endl;
    energy = kineticEnergy + potentialEnergy;
}

void StatisticsCalculator::calculateTemperature(){
    temperature = 2./(3.*mySystem->unFrozen)*kineticEnergy;
}

void StatisticsCalculator::calculatePressureFast(){
    //pressure = mySystem->fastPressure;
    pressure = mySystem->unFrozen/volume*temperature + (1./(3.*volume))*pressure;
}

void StatisticsCalculator::calculatePressure(){

}

void StatisticsCalculator::calculateDisplacement(){
    displacement = 0;
    for(int i=0; i<mySystem->totalAtoms; i++){
        displacement += dot(mySystem->atomList[i]->getRealPos() - mySystem->atomList[i]->getInitialPos(), mySystem->atomList[i]->getRealPos() - mySystem->atomList[i]->getInitialPos());
    }
    displacement /= mySystem->unFrozen;
}

void StatisticsCalculator::calculateVolume(){

}

void StatisticsCalculator::sampleStats(){
    calculateEnergy();
    calculateTemperature();
    calculatePressureFast();
    calculateDisplacement();
    calculateVolume();
}
