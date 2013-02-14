#include "system.h"
using namespace std;
using namespace arma;
System::System()
{
    nx = 8;
    ny = 8;
    nz = 8;
    b = 1.545;
    atomsPerGridPoint = 4;
    totalAtoms = atomsPerGridPoint*nx*ny*nz;
    for(int i=0; i<totalAtoms;i++){
        atomList.push_back(new Atom());
        atomList[i]->setSystemIndex(i);
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
    cellSize = 2.0;
    cellsInXDir = (int) floor(nx*b/cellSize);
    cellsInYDir = (int) floor(ny*b/cellSize);
    cellsInZDir = (int) floor(nz*b/cellSize);
    totalCells = cellsInXDir*cellsInYDir*cellsInZDir;
    cellSize = nx*b/cellsInXDir; //updata cell size to actual
    cout << cellSize << " " << totalCells << endl;
    vec3 posvec;
    int cellNumber = 0;
    int neighbour = 0;
    int x;
    int y;
    int z;
    cout << totalCells << endl;
    for(int k=0;k<cellsInZDir;k++){
        for(int j=0;j<cellsInYDir;j++){
            for(int i=0;i<cellsInXDir;i++){
                posvec << i*cellSize << j*cellSize << k*cellSize;
                cellNumber = k*cellsInXDir*cellsInYDir + j*cellsInXDir + i;
                cellList.push_back(new Cell(posvec, cellSize, cellNumber));
                //create neighbourcells
                x = 1;
                for(y=-1;y<2;y++){
                    for(z=-1;z<2;z++){
                        neighbour = ((k+z+cellsInZDir) % cellsInZDir)*cellsInXDir*cellsInYDir + ((j+y+cellsInYDir) % cellsInYDir)*cellsInXDir + ((i+x+cellsInXDir) % cellsInXDir);
                        if(neighbour != cellNumber){
                            cellList[cellNumber]->neighbourList.push_back(neighbour);
                        }
                    }
                }
                x = 0;
                y = 1;
                for(z=-1;z<2;z++){
                    neighbour = ((k+z+cellsInZDir) % cellsInZDir)*cellsInXDir*cellsInYDir + ((j+y+cellsInYDir) % cellsInYDir)*cellsInXDir + ((i+x+cellsInXDir) % cellsInXDir);
                    if(neighbour != cellNumber){
                        cellList[cellNumber]->neighbourList.push_back(neighbour);
                    }
                }
                y = 0;
                z = 1;
                neighbour = ((k+z+cellsInZDir) % cellsInZDir)*cellsInXDir*cellsInYDir + ((j+y+cellsInYDir) % cellsInYDir)*cellsInXDir + ((i+x+cellsInXDir) % cellsInXDir);
                if(neighbour != cellNumber){
                    cellList[cellNumber]->neighbourList.push_back(neighbour);
                }
            }
        }
    }

    noOfTimeSteps = 10;
    setPosFCC();
    setVelNormal();
    cout << "placing atoms in initial cells... ";
    placeAtomsInCells();
    cout << "done!" << endl;
    cout << "Calculating first time forces... ";
    calculateForcesCellsLJMIC();
    cout << "done!" << endl;


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
                  vel(0) << " " << vel(1) << " " << vel(2) << " " << cellList[atomList[i]->getCellNumber()]->colorIndex << " " <<std::endl;

    }
    myfile.close();
}


void System::setVelNormal(){
    cout << "Setting normally distributed velocities... ";
    vec3 randomvec;
    vec3 totaldrift = zeros(3);
    for(int i=0; i<totalAtoms;i++){
        randomvec = randn<vec>(3);
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
    placeAtomsInCells();
    calculateForcesCellsLJMIC();
    for(int i=0; i<totalAtoms;i++){
        atomList[i]->setVel(atomList[i]->getVel() + forces.col(i)*dt/(2*mass));
    }
}

void System::calculateForcesNull(){
    vec3 tempvec = zeros(3);
    for(int i=0; i<totalAtoms;i++){
        forces.col(i) = tempvec;
    }
}

void System::calculateForcesLJMIC(){
    // leonard-jones force using minimal image convention
    calculateForcesNull(); //zero out the forces vector
    vec3 singlePair;
    for(int i=0; i<totalAtoms;i++){
        for(int j=i+1;j<totalAtoms;j++){
            singlePair = (atomList[i]->getPos() - atomList[j]->getPos());
            singlePairForces(singlePair);
            forces.col(i) += singlePair;
            forces.col(j) -= singlePair;
        }
    }
}

void System::singlePairForces(vec3& singlePair){
    double r2;
    double r6;
    double r8;
    double factor;
    for (int j=0; j<3; j++){
        singlePair(j) = (abs(singlePair(j)) < abs(singlePair(j) + nx*b))*(abs(singlePair(j)) < abs(singlePair(j) - nx*b))*singlePair(j)
                + (abs(singlePair(j) + nx*b) < abs(singlePair(j)))*(abs(singlePair(j) + nx*b) < abs(singlePair(j) - nx*b))*(singlePair(j)+nx*b)
                + (abs(singlePair(j) - nx*b) < abs(singlePair(j)))*(abs(singlePair(j) - nx*b) < abs(singlePair(j) + nx*b))*(singlePair(j)-nx*b);
    }
    r2 = max(singlePair(0)*singlePair(0) +  singlePair(1)*singlePair(1) + singlePair(2)*singlePair(2),0.81);
    //r2 = max(dot(singlePair,singlePair),0.81);
    r6 = r2*r2*r2;
    r8 = r2*r6;
    factor = (24./r8)*(2/r6 - 1);
    singlePair(0) = factor*singlePair(0);
    singlePair(1) = factor*singlePair(1);
    singlePair(2) = factor*singlePair(2);
    //singlePair = singlePair*(24./r8)*(2/r6 - 1);
}

void System::runSimulation(){
    ostringstream number;
    string filename = "results";
    string saveFilename;
    string extension = ".xyz";
    vmdPrintSystem(saveFilename);
    for(int i=0; i<noOfTimeSteps;i++){
        cout << "solving for time step: " << i << " out of: " << noOfTimeSteps << "...";

        //some jargon to make file extensions 000, 001, 002, etc.
        filename = "results";
        number << setw(3) << setfill('0') << i;
        saveFilename = number.str();
        number.str(""); //TODO: This is bad use of ostringstream, not meant to be reused!
        saveFilename = filename.append(saveFilename).append(extension);
        //jargon done!
        vmdPrintSystem(saveFilename);
        timeEvolve();

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
        xCell = (int) (posvec(0)/cellSize);
        yCell = (int) (posvec(1)/cellSize);
        zCell = (int) (posvec(2)/cellSize);
        oldCellNumber = atomList[i]->getCellNumber();

        cellNumber = zCell*cellsInYDir*cellsInXDir + yCell*cellsInXDir + xCell;
        if(cellNumber != oldCellNumber){
            if(oldCellNumber != -1){
                cellList[oldCellNumber]->atomsInCell.remove(atomList[i]);
            }
        cellList[cellNumber]->atomsInCell.push_back(atomList[i]);
        atomList[i]->setCellNumber(cellNumber);
        }
    }
}


void System::calculateForcesCellsLJMIC(){
    // leonard-jones force using minimal image convention
    calculateForcesNull(); //zero out the forces vector
    vec3 singlePair;
    int i;
    int thisIndex;
    int otherIndex;
    list<Atom*>::iterator dummyk; //IM NOT EVEN JOKING WTF IS THIS
    int jint = 0; //SRSLY THIS IS WHAT I HAVE TO DEAL WITH
    list<Atom*>::iterator j;
    list<Atom*>::iterator k;
    int neighbourPos;
    for(int zCellPos=0; zCellPos<cellsInXDir;zCellPos++){
        for(int yCellPos=0; yCellPos<cellsInYDir;yCellPos++){
            for(int xCellPos=0; xCellPos<cellsInYDir;xCellPos++){
                i = zCellPos*cellsInXDir*cellsInYDir + yCellPos*cellsInXDir + xCellPos;
                for(j=cellList[i]->atomsInCell.begin(); j!= cellList[i]->atomsInCell.end(); ++j){
                    //particles in own cell
                    thisIndex = (*j)->getSystemIndex();
                    dummyk = cellList[i]->atomsInCell.begin();
                    advance(dummyk, jint+1);
                    for(k=dummyk; k != cellList[i]->atomsInCell.end(); ++k){
                        otherIndex = (*k)->getSystemIndex();
                        singlePair = (atomList[thisIndex]->getPos() - atomList[otherIndex]->getPos());
                        singlePairForces(singlePair);
                        forces.col(thisIndex) += singlePair;
                        forces.col(otherIndex) -= singlePair;
                    }
                    jint ++;

                    //particles in neighbouring cells
                    //each cell points to 13 neighbouring cells, to avoid double counting
                    //the cell points to cells which have (x,y,z)-indices increased by 1 (mod N)
                    //each cell is initialized with a neighbour vector<int> called neighbourList
                    for(int neighbour = 0; neighbour<cellList[i]->neighbourList.size(); neighbour++){
                        neighbourPos = cellList[i]->neighbourList[neighbour];
                        for(k=cellList[neighbourPos]->atomsInCell.begin(); k!= cellList[neighbourPos]->atomsInCell.end(); ++k){
                            otherIndex = (*k)->getSystemIndex();
                            singlePair = (atomList[thisIndex]->getPos() - atomList[otherIndex]->getPos());
                            singlePairForces(singlePair);
                            forces.col(thisIndex) += singlePair;
                            forces.col(otherIndex) -= singlePair;
                        }
                    }
                }
            }
        }
    }
}
