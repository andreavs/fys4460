#include <iostream>
#include <stdlib.h>
#include "particle.h"
#include "cell.h"
#include <string>
#include <iomanip>
#include <fstream>

using namespace std;


int main()
{
    cout << "Hello World!" << endl;

    int Nx = 8; int Ny = 8; int Nz = 8;
    Cell *argongrid = new Cell(Nx,Ny,Nz);
    double b = 5.260; // Aangstroms
    int counter = 0;

    for(int i=0;i<Nx;i++){
        for(int j=0;j<Ny;j++){
            for(int k=0;k<Nz;k++){
                argongrid->cellContents[4*(Ny*Nz*i + Nz*j + k)] = new Particle(1+b*i,b*j,b*k,0.0,0.0,0.0);
                argongrid->cellContents[4*(Ny*Nz*i + Nz*j + k)]->setMoleculeType("Ar");

                argongrid->cellContents[4*(Ny*Nz*i + Nz*j + k)+1] = new Particle(b*(i+0.5),b*(j+0.5),b*k,0.0,0.0,0.0);
                argongrid->cellContents[4*(Ny*Nz*i + Nz*j + k)+1]->setMoleculeType("Ar");

                argongrid->cellContents[4*(Ny*Nz*i + Nz*j + k)+2] = new Particle(b*i,b*(j+0.5),b*(k+0.5),0.0,0.0,0.0);
                argongrid->cellContents[4*(Ny*Nz*i + Nz*j + k)+2]->setMoleculeType("Ar");

                argongrid->cellContents[4*(Ny*Nz*i + Nz*j + k)+3] = new Particle(b*(i+0.5),b*j,b*(k+0.5),0.0,0.0,0.0);
                argongrid->cellContents[4*(Ny*Nz*i + Nz*j + k)+3]->setMoleculeType("Ar");
            }
        }
    }
    string a = "hello";
    cout << a.append(" cool") << endl;
    std::string filename = "results.xyz";
    argongrid->vmdPrintCell(filename);
    cout << filename << endl;


    return 0;



}

