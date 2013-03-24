#include "atom.h"
#include "configreader.h"
using namespace arma;
using namespace std;


Atom::Atom(string an, bool freeze)
{
    atomName = an;
    pos = zeros<vec>(3);
    vel = zeros<vec>(3);
    cellNumber = -1;
    mass = 1.0;
    isFrozen = freeze;
}

Atom::Atom(ConfigReader *cfgReader, string an, bool freeze){
    atomName = an;
    pos = zeros<vec>(3);
    vel = zeros<vec>(3);
    cellNumber = -1;
    isFrozen = freeze;

}
