#include "system.h"
using namespace std;
using namespace arma;
System::System()
{
    nx = 12;
    ny = 12;
    nz = 12;
    b = 1.545;
    atomsPerGridPoint = 4;
    totalAtoms = atomsPerGridPoint*nx*ny*nz;
    for(int i=0; i<totalAtoms;i++){
        atomList.push_back(new Atom());
    }
    forces = zeros(3,totalAtoms);
    mass = 1.0; // mass in argon mass units
    time = 0.0;
    dt = 0.01; // time step
    F0 = 1.0; // force thing in md units
    E0 = 1.0; // energy
    T0 = 1.0; // temperature
    kb = E0/T0;
    sigma = sqrt(kb*T0/mass);
    cellSize = 3*sigma;
    cellsInXDir = (int) floor(nx*b/cellSize+1);
    cellsInYDir = (int) floor(ny*b/cellSize+1);
    cellsInZDir = (int) floor(nz*b/cellSize+1);
    totalCells = cellsInXDir*cellsInYDir*cellsInZDir;
    vec3 posvec;
    for(int i=0;i<cellsInXDir;i++){
        for(int j=0;j<cellsInYDir;j++){
            for(int k=0;k<cellsInZDir;k++){
                posvec << i*cellSize << j*cellSize << k*cellSize;
                cellList.push_back(new Cell(posvec, cellSize));
            }
        }
    }
    noOfTimeSteps = 500;
    setPosFCC();
    setVelNormal();

    placeAtomsInCells();
    calculateForcesLJMIC();

}


void System::vmdPrintSystem(std::string filename)
{
    std::ofstream myfile;
    myfile.open(filename.c_str());
    vec3 pos;
    vec3 vel;
    std::string name;
    myfile << totalAtoms << std::endl;
    myfile << "This line has not unintentionally been left unblank" << std::endl;
    for(int i=0;i<totalAtoms;i++){
        name = atomList[i]->getName();
        pos = atomList[i]->getPos();
        vel = atomList[i]->getVel();

        myfile << name << " " << pos(0) << " " << pos(1) << " " << pos(2) << " " <<
                  vel(0) << " " << vel(1) << " " << vel(2) << " " << std::endl;

    }
    myfile.close();
}


void System::setVelNormal(){
    cout << "Setting normally distributed velocities... ";
    vec3 randomvec;
    vec3 totaldrift = zeros(3);
    for(int i=0; i<totalAtoms;i++){
        randomvec = sigma*randn<vec>(3);
        atomList[i]->setVel(randomvec);
        totaldrift += randomvec;
    }
    totaldrift = totaldrift/totalAtoms;
    for(int i=0; i<totalAtoms;i++){
        atomList[i]->setVel(atomList[i]->getVel() - totaldrift);
    }
    cout << "done" << endl;
}


void System::setPosFCC(){
    cout << "Configuring face centered cube grid... ";
    vec3 posvec;
    for(int i=0;i<nx;i++){
        for(int j=0;j<ny;j++){
            for(int k=0;k<nz;k++){
                posvec << b*i << b*j << b*k;
                atomList[4*(ny*nz*k + nz*j + i)]->setPos(posvec);

                posvec << b*(i+0.5) << b*(j+0.5) << b*k;
                atomList[4*(ny*nz*k + nz*j + i)+1]->setPos(posvec);

                posvec << b*i << b*(j+0.5) << b*(k+0.5);
                atomList[4*(ny*nz*k + nz*j + i)+2]->setPos(posvec);

                posvec << b*(i+0.5) << b*j << b*(k+0.5);
                atomList[4*(ny*nz*k + nz*j + i)+3]->setPos(posvec);

            }
        }
    }
    cout << "done" << endl;
}


void System::timeEvolve(){

    for(int i=0; i<totalAtoms;i++){
        forces.col(i) = atomList[i]->getVel() + forces.col(i)*dt/(2*mass);
        atomList[i]->setPos(atomList[i]->getPos() + forces.col(i)*dt);
    }
    time += dt;
    periodicBoundaries();
    calculateForcesLJMIC();
    for(int i=0; i<totalAtoms;i++){
        atomList[i]->setVel(atomList[i]->getVel() + forces.col(i)*dt/(2*mass));
    }
    placeAtomsInCells();
}

void System::calculateForcesNull(){
    vec3 tempvec = zeros(3);
    for(int i=0; i<totalAtoms;i++){
        forces.col(i) = tempvec;
    }
}

void System::calculateForcesLJMIC(){
    // leonard-jones force using minimal image convention
    vec3 singlePair;
    double r;
    for(int i=0; i<totalAtoms;i++){
        for(int j=i+1;j<totalAtoms;j++){
            singlePair = (atomList[i]->getPos() - atomList[j]->getPos());
            for (int j=0; j<3; j++){
                singlePair(j) = (abs(singlePair(j)) < abs(singlePair(j) + nx*b))*(abs(singlePair(j)) < abs(singlePair(j) - nx*b))*singlePair(j)
                        + (abs(singlePair(j) + nx*b) < abs(singlePair(j)))*(abs(singlePair(j) + nx*b) < abs(singlePair(j) - nx*b))*(singlePair(j)+nx*b)
                        + (abs(singlePair(j) - nx*b) < abs(singlePair(j)))*(abs(singlePair(j) - nx*b) < abs(singlePair(j) + nx*b))*(singlePair(j)-nx*b);
            }
            //cout << atomList[j]->getPos()r << endl;
            //singlePair.print();
            r = max(norm(singlePair,2),0.8);
            //cout << r << endl;
            singlePair = singlePair*(24./pow(r,8))*(2/pow(r,6) - 1);
            forces.col(i) += singlePair;
            forces.col(j) -= singlePair;
        }
    }

}

void System::runSimulation(){
    ostringstream number;
    string filename = "results";
    string saveFilename;
    string extension = ".xyz";
    vmdPrintSystem(saveFilename);
    for(int i=0; i<noOfTimeSteps;i++){
        cout << "solving for time step: " << i << " out of: " << noOfTimeSteps << "...";
        timeEvolve();

        //some jargon to make file extensions 000, 001, 002, etc.
        filename = "results";
        number << setw(3) << setfill('0') << i;
        saveFilename = number.str();
        number.str(""); //TODO: This is bad use of ostringstream, not meant to be reused!
        saveFilename = filename.append(saveFilename).append(extension);
        //jargon done!

        vmdPrintSystem(saveFilename);
        cout << " done!" << endl;
    }
}

void System::periodicBoundaries(){
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

void System::placeAtomsInCells(){
    int xCell;
    int yCell;
    int zCell;
    int cellNumber;
    int oldCellNumber;
    vec posvec;
    for(int i=0;i<totalAtoms;i++){
        posvec = atomList[i]->getPos();
        xCell = (int) posvec(0)/cellSize;
        yCell = (int) posvec(1)/cellSize;
        zCell = (int) posvec(2)/cellSize;
        oldCellNumber = atomList[i]->getCellNumber();
        cellNumber = zCell*cellsInYDir*cellsInXDir + yCell*cellsInXDir + xCell;
        if(cellNumber != oldCellNumber){
            cout << i << endl;
            if(oldCellNumber != -1){
                cellList[oldCellNumber]->atomsInCell.erase(cellList[oldCellNumber]->atomsInCell.begin()
                                                           + atomList[i]->getPosInCell()-1);
            }
            cout << i << endl;
            cout << oldCellNumber << " " << cellNumber << " " << cellList.size() << endl;
            cellList[cellNumber]->atomsInCell.push_back(atomList[i]);
            atomList[i]->setCellNumber(cellNumber);
            atomList[i]->setPosInCell(cellList[cellNumber]->atomsInCell.size()-1);
        }
    }
}


void System::calculateForcesCellsLJMIC(){

}
