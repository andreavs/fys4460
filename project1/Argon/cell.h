#ifndef CELL_H
#define CELL_H
#include <iostream>
#include <string>
#include <stdio.h>
#include "particle.h"
#include <iomanip>
#include <fstream>



class Cell
{
public:
    Cell(int lengthx, int lengthy, int lengthz);

    void vmdPrintCell(std::string filename);
    Particle **cellContents;

    int* getSize();

private:
    int Nx; int Ny; int Nz;
};

#endif // CELL_H
