#include "cell.h"

Cell::Cell(int lengthx,int lengthy,int lengthz)
{
    Nx = lengthx;
    Ny = lengthy;
    Nz = lengthz;
    //std::cout << Ny << Nz << std::endl;
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
    //double *pos = new double[3];
    double *vel = new double[3];
    std::string name;
    int prod = 4*Nx * Ny * Nz;
    myfile << prod << std::endl;
    myfile << "This line has not unintentionally been left unblank" << std::endl;
    std::cout << prod << std::endl;


    //std::cout << "balle43234"  << std::endl;
    //double*pos = Cell::cellContents[0]->getPos(pos);
    //std::cout << "balle43235" << " "<< pos[1] << std::endl;


    for(int i=0;i<prod;i++){
        //std::cout << "balle " << pos[1] << std::endl;
        name = cellContents[i]->getMoleculeType();
        std::cout << cellContents[i]->pos[0];
        //std::cout << "balle43234" << pos[1] << std::endl;
        //pos = Cell::cellContents[i]->getPos(pos);
        std::cout << "balle" << i << std::endl;
        //vel = Cell::cellContents[i]->getVel();
        //std::cout << "balle" << pos[1] << std::endl;
        myfile << name << " " << cellContents[i]->pos[0] << " " << cellContents[i]->pos[1] << " " << cellContents[i]->pos[2] << " " <<
                  cellContents[i]->vel[0] << " " << cellContents[i]->vel[1] << " " << cellContents[i]->vel[2] << " " << std::endl;

    }
    myfile.close();
}
