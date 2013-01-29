#include "cell.h"

Cell::Cell(int lengthx,int lengthy,int lengthz)
{
    Nx = lengthx;
    Ny = lengthy;
    Nz = lengthz;
    std::cout << Ny << Nz << std::endl;
    Cell::cellContents = new Particle*[Nx*Ny*Nz*4];

}

int* Cell::getSize()
{
    int* size = new int[3];
    size[0] = Nx; size[1] = Ny; size[2] = Nz;
    return size;
}


/////////////////////////////////
// vmdPrintcell(filename)      //
// prints the content of the   //
// cell in .xyz format.        //
/////////////////////////////////

void Cell::vmdPrintCell(std::string filename)
{
    std::ofstream myfile;
    //filename = filename.append(".xyz");
    myfile.open(filename.c_str());
    double *pos = new double[3];
    double *vel = new double[3];
    std::string name;
    int prod = 4*Nx * Ny * Nz;
    myfile << prod << std::endl;
    myfile << "This line has not unintentionally been left unblank" << std::endl;
    std::cout << prod << std::endl;
    for(int i=0;i<prod;i++){
        std::cout << "balle" << pos[1] << std::endl;
        name = Cell::cellContents[i]->getMoleculeType();
        std::cout << "balle" << pos[1] << std::endl;
        //pos = Cell::cellContents[i]->getPos();
        std::cout << "balle" << pos[1] << std::endl;
        //vel = Cell::cellContents[i]->getVel();
        std::cout << "balle" << pos[1] << std::endl;
        myfile << name << " " << pos[0] << " " << pos[1] << " " << pos[2] << " " <<
                  vel[0] << " " << vel[1] << " " << vel[2] << " " << std::endl;

    }
    myfile.close();
}
