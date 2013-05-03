#include "simulation.h"
#include "configreader.h"
#include "system.h"
#include "atom.h"
#include "cell.h"
#include "statisticscalculator.h"
#include "tools.cpp"



using namespace std;
using namespace arma;

System::System(ConfigReader *cfgReader)
{

    try
    {
      //fcc spesifics:
      nx = cfgReader->cfg.lookup("nx");
      ny = cfgReader->cfg.lookup("ny");
      nz = cfgReader->cfg.lookup("nz");
      atomsPerGridPoint = cfgReader->cfg.lookup("atomsPerGridPoint");
      createFCC = cfgReader->cfg.lookup("createFCC");
      setVelNormalbool = cfgReader->cfg.lookup("setVelNormal");
      b = cfgReader->cfg.lookup("b");

      //read initial from .xyz file instead:
      readInitialFromFile = cfgReader->cfg.lookup("readInitialFromFile");
      string exp = cfgReader->cfg.lookup("experiment");
      experiment = exp;
      //cout << experiment << endl;
      // other info
      cellSize = cfgReader->cfg.lookup("cellSize");
      time = cfgReader->cfg.lookup("time");
      tempInKelvin = cfgReader->cfg.lookup("tempInKelvin");
      //string an = cfgReader->cfg.lookup("particleName");
      //cout << an << endl;

    }
    catch(const SettingNotFoundException &nfex)
    {
      cerr << "No 'name' setting in configuration file when read from System." << endl;
    }
    temperature = tempInKelvin/119.74;
    stddev = sqrt(temperature);
    experiment = "experiments/" + experiment + "/results/lastState.xyz";


    if(createFCC){
        totalAtoms = atomsPerGridPoint*nx*ny*nz;
        for(int i=0; i<totalAtoms;i++){
            atomList.push_back(new Atom());
            atomList[i]->setSystemIndex(i);
        }
        setPosFCC();
        if(setVelNormalbool){
            setVelNormal();
        }
        else if(setVelUniformbool){
            setVelUniform();
        }
    }
    else if(readInitialFromFile){
        ifstream myfile (experiment);
        string line;
        vector<string> contents;
        cout << experiment << endl;
        string name;
        vec3 pos = zeros(3);
        vec3 vel = zeros(3);
        bool freeze;
        cout << "starting from experiment " << experiment << endl;
        if (myfile.is_open())
        {
            getline(myfile,line);
            totalAtoms = atoi(line.c_str());
            getline(myfile,line);
            cout << totalAtoms << endl;

            int counter = 0;
            for(int i=0;i<totalAtoms;i++)
            {
                getline (myfile,line);
                contents = split(line, ' ', contents);
                name = contents[0];
                pos(0) = atof(contents[1].c_str());
                pos(1) = atof(contents[2].c_str());
                pos(2) = atof(contents[3].c_str());
                vel(0) = atof(contents[4].c_str());
                vel(1) = atof(contents[5].c_str());
                vel(2) = atof(contents[6].c_str());
                freeze = to_bool(contents[10]);
                atomList.push_back(new Atom(name, freeze));
                //cout << counter << endl;
                atomList[i]->setPos(pos);
                atomList[i]->initialPos = pos;
                atomList[i]->realPos = pos;
                atomList[i]->setVel(vel);
                atomList[i]->setSystemIndex(i);
                contents.clear();

                //cout << i << endl;
                counter++;
            }
        //cout << counter << endl;
        //totalAtoms = atomList.size();
        //cout << totalAtoms << endl;
        myfile.close();
        }
        else cout << "Unable to open file";
       }


    cout << "setting up computation cells... ";
    setupCells();
    cout << "done!" << endl;

    unFrozen = 0;
    for(int i=0;i<totalAtoms;i++){
        if(!atomList[i]->isFrozen){
            unFrozen++;
        }
    }


    placeAtomsInCells();
}

void System::setupCells(){
    cellsInXDir = (int) floor(nx*b/cellSize);
    cellsInYDir = (int) floor(ny*b/cellSize);
    cellsInZDir = (int) floor(nz*b/cellSize);
    totalCells = cellsInXDir*cellsInYDir*cellsInZDir;
    cellSize = nx*b/cellsInXDir; //updata cell size to actual
    cout << "real cellsize: " << cellSize << " total cells: " << totalCells << "...";
    vec3 posvec;
    int cellNumber = 0;
    int neighbour = 0;
    volume = (b*nx)*(b*ny)*(b*nz);
    int x;
    int y;
    int z;
    vec3 minimalImage;
    // lots of farblegarble to make neighbourlists incoming!
    for(int k=0;k<cellsInZDir;k++){
        for(int j=0;j<cellsInYDir;j++){
            for(int i=0;i<cellsInXDir;i++){
                posvec << i*cellSize << j*cellSize << k*cellSize;
                cellNumber = k*cellsInXDir*cellsInYDir + j*cellsInXDir + i;
                cellList.push_back(new Cell(posvec, cellSize, cellNumber));
            }
        }
    }
    //create neighbourcells
    int n [3];
    n[0] = nx; n[1] = ny; n[2] = nz;
    for(int k=0;k<cellsInZDir;k++){
        for(int j=0;j<cellsInYDir;j++){
            for(int i=0;i<cellsInXDir;i++){
                cellNumber = k*cellsInXDir*cellsInYDir + j*cellsInXDir + i;
                x = 1;
                for(y=-1;y<2;y++){
                    for(z=-1;z<2;z++){
                        neighbour = ((k+z+cellsInZDir) % cellsInZDir)*cellsInXDir*cellsInYDir + ((j+y+cellsInYDir) % cellsInYDir)*cellsInXDir + ((i+x+cellsInXDir) % cellsInXDir);
                        if(neighbour != cellNumber){
                            cellList[cellNumber]->neighbourList.push_back(neighbour);
                            minimalImage = cellList[neighbour]->cellPos - cellList[cellNumber]->cellPos;
                            for (int j=0; j<3; j++){
                                minimalImage(j) = (abs(minimalImage(j)) < abs(minimalImage(j) + n[j]*b))*(abs(minimalImage(j)) < abs(minimalImage(j) - n[j]*b))*minimalImage(j)
                                        + (abs(minimalImage(j) + n[j]*b) < abs(minimalImage(j)))*(abs(minimalImage(j) + n[j]*b) < abs(minimalImage(j) - n[j]*b))*(minimalImage(j)+n[j]*b)
                                        + (abs(minimalImage(j) - n[j]*b) < abs(minimalImage(j)))*(abs(minimalImage(j) - n[j]*b) < abs(minimalImage(j) + n[j]*b))*(minimalImage(j)-n[j]*b)
                                        - minimalImage(j);
                            }
                            cellList[cellNumber]->neighbourVectors.push_back(minimalImage);
                        }
                    }
                }
                x = 0;
                y = 1;
                for(z=-1;z<2;z++){
                    neighbour = ((k+z+cellsInZDir) % cellsInZDir)*cellsInXDir*cellsInYDir + ((j+y+cellsInYDir) % cellsInYDir)*cellsInXDir + ((i+x+cellsInXDir) % cellsInXDir);
                    if(neighbour != cellNumber){
                        cellList[cellNumber]->neighbourList.push_back(neighbour);
                        minimalImage = cellList[neighbour]->cellPos - cellList[cellNumber]->cellPos;
                        for (int j=0; j<3; j++){
                            minimalImage(j) = (abs(minimalImage(j)) < abs(minimalImage(j) + n[j]*b))*(abs(minimalImage(j)) < abs(minimalImage(j) - n[j]*b))*minimalImage(j)
                                    + (abs(minimalImage(j) + n[j]*b) < abs(minimalImage(j)))*(abs(minimalImage(j) + n[j]*b) < abs(minimalImage(j) - n[j]*b))*(minimalImage(j)+n[j]*b)
                                    + (abs(minimalImage(j) - n[j]*b) < abs(minimalImage(j)))*(abs(minimalImage(j) - n[j]*b) < abs(minimalImage(j) + n[j]*b))*(minimalImage(j)-n[j]*b)
                                    - minimalImage(j);
                        }
                        cellList[cellNumber]->neighbourVectors.push_back(minimalImage);
                    }
                }
                y = 0;
                z = 1;
                neighbour = ((k+z+cellsInZDir) % cellsInZDir)*cellsInXDir*cellsInYDir + ((j+y+cellsInYDir) % cellsInYDir)*cellsInXDir + ((i+x+cellsInXDir) % cellsInXDir);
                if(neighbour != cellNumber){
                    cellList[cellNumber]->neighbourList.push_back(neighbour);
                    minimalImage = cellList[neighbour]->cellPos - cellList[cellNumber]->cellPos;
                    for (int j=0; j<3; j++){
                        minimalImage(j) = (abs(minimalImage(j)) < abs(minimalImage(j) + n[j]*b))*(abs(minimalImage(j)) < abs(minimalImage(j) - n[j]*b))*minimalImage(j)
                                + (abs(minimalImage(j) + n[j]*b) < abs(minimalImage(j)))*(abs(minimalImage(j) + n[j]*b) < abs(minimalImage(j) - n[j]*b))*(minimalImage(j)+n[j]*b)
                                + (abs(minimalImage(j) - n[j]*b) < abs(minimalImage(j)))*(abs(minimalImage(j) - n[j]*b) < abs(minimalImage(j) + n[j]*b))*(minimalImage(j)-n[j]*b)
                                - minimalImage(j);
                    }
                    cellList[cellNumber]->neighbourVectors.push_back(minimalImage);
                }
            }
        }
    }
    // farblegarble done!
}

void System::setVelNormal(){
    cout << "Setting normally distributed velocities... ";
    vec3 randomvec;
    vec3 totaldrift = zeros(3);
    for(int i=0; i<totalAtoms;i++){
        randomvec = stddev*randn<vec>(3);
        atomList[i]->setVel(randomvec);
        totaldrift += randomvec;
    }
    totaldrift = totaldrift/totalAtoms;
    for(int i=0; i<totalAtoms;i++){
        atomList[i]->setVel(atomList[i]->getVel() - totaldrift);
    }
    cout << "done" << endl;
}

void System::setVelUniform(){
    cout << "Setting uniformly distributed velocities... ";
    vec3 randomvec;
    vec3 totaldrift = zeros(3);
    for(int i=0; i<totalAtoms;i++){
        randomvec = 2*stddev*(randu<vec>(3)-0.5);
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
                atomList[4*(ny*nz*k + nz*j + i)]->setInitialPos(posvec);
                atomList[4*(ny*nz*k + nz*j + i)]->setRealPos(posvec);

                posvec << b*(i+0.5) << b*(j+0.5) << b*k;
                atomList[4*(ny*nz*k + nz*j + i)+1]->setPos(posvec);
                atomList[4*(ny*nz*k + nz*j + i)+1]->setInitialPos(posvec);
                atomList[4*(ny*nz*k + nz*j + i)+1]->setRealPos(posvec);

                posvec << b*i << b*(j+0.5) << b*(k+0.5);
                atomList[4*(ny*nz*k + nz*j + i)+2]->setPos(posvec);
                atomList[4*(ny*nz*k + nz*j + i)+2]->setInitialPos(posvec);
                atomList[4*(ny*nz*k + nz*j + i)+2]->setRealPos(posvec);

                posvec << b*(i+0.5) << b*j << b*(k+0.5);
                atomList[4*(ny*nz*k + nz*j + i)+3]->setPos(posvec);
                atomList[4*(ny*nz*k + nz*j + i)+3]->setInitialPos(posvec);
                atomList[4*(ny*nz*k + nz*j + i)+3]->setRealPos(posvec);

            }
        }
    }
    cout << "done" << endl;
}

void System::placeAtomsInCells(){
    int xCell;
    int yCell;
    int zCell;
    int cellNumber;
    int oldCellNumber;
    //int cells = cellsInYDir*cellsInXDir*cellsInZDir;

    vec posvec;
    for(int i=0;i<totalAtoms;i++){
        //cout << i << endl;
        posvec = atomList[i]->getPos();

        xCell = (int) (posvec(0)/cellSize);
        yCell = (int) (posvec(1)/cellSize);
        zCell = (int) (posvec(2)/cellSize);
        oldCellNumber = atomList[i]->getCellNumber();
        cellNumber = zCell*cellsInYDir*cellsInXDir + yCell*cellsInXDir + xCell;
        //cout << cellNumber << endl;
        if(cellNumber != oldCellNumber){
            if(oldCellNumber != -1){

                cellList[oldCellNumber]->atomsInCell.remove(atomList[i]);
            }
        cellList[cellNumber]->atomsInCell.push_back(atomList[i]);
        atomList[i]->setCellNumber(cellNumber);
        }
    }
}

void System::vmdPrintSystem(std::string filename)
{
    std::ofstream myfile;
    myfile.open(filename.c_str());
    vec3 pos;
    vec3 vel;
    vec3 displacementvec;
    double displacement;
    std::string name;
    myfile << totalAtoms << std::endl;
    myfile << "This line has not unintentionally been left unblank" << std::endl;
    for(int i=0;i<totalAtoms;i++){
        name = atomList[i]->getName();
        pos = atomList[i]->getPos();
        vel = atomList[i]->getVel();
        displacementvec = atomList[i]->getRealPos() - atomList[i]->getInitialPos();

        displacement = dot(displacementvec,displacementvec);
        myfile << name << " " << pos(0) << " " << pos(1) << " " << pos(2) << " " <<
                  vel(0) << " " << vel(1) << " " << vel(2) << " " << cellList[atomList[i]->getCellNumber()]->colorIndex << " "
               << displacement << " " << norm(vel,2) << " " << atomList[i]->isFrozen << " " << atomList[i]->atomPressure << " " << std::endl;
    }
    myfile.close();
}
