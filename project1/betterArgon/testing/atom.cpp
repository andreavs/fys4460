#include "atom.h"
using namespace arma;
using namespace std;


Atom::Atom()
{
    atomName = "Ar";
    pos = zeros<vec>(3);
    vel = zeros<vec>(3);
    cellNumber = -1; //initialize with impossible cell to make first test
    posInCell = -1;


}
