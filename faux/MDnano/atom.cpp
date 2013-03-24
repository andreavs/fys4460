#include "atom.h"
#include "configreader.h"
using namespace arma;
using namespace std;


Atom::Atom()
{
    atomName = "Ar";
    pos = zeros<vec>(3);
    vel = zeros<vec>(3);
    cellNumber = -1;
    mass = 1.0;
}

Atom::Atom(ConfigReader *cfgReader){
    try
    {
        atomName = cfgReader->cfg.lookupValue("particleName", atomName);

    }
    catch(const SettingNotFoundException &nfex)
    {
      cerr << "No 'name' setting in configuration file when read from Atom." << endl;
    }
    pos = zeros<vec>(3);
    vel = zeros<vec>(3);
    cellNumber = -1;

}
