#include "atom.h"
using namespace arma;
using namespace std;


Atom::Atom()
{
    atomName = "Ar";
    pos = zeros<vec>(3);
    vel = zeros<vec>(3);


}
