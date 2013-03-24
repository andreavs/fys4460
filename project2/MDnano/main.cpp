
#include <iostream>
#include <string>
#include <iostream>
#include <iomanip>
#include <cstdlib>

#include <simulation.h>
#include <cell.h>
#include <atom.h>
#include <system.h>


#include<libconfig.h++>


using namespace std;
using namespace libconfig;

int main(int argc, char* argv[])
{
    cout << "Hello World!" << endl;
    string fn = argv[1];
    string filename;
    try{
        filename = "experiments/" + fn + "/" + fn + ".cfg";
    }
    catch(int e){
        cout << "file not found" << endl;
        return 0;
    }
    Simulation *mySim = new Simulation(fn);
    mySim->runSimulation();

//    Config cfg;

//    char *fn = new char[1000];
//    fn = "test.cfg";
//    cout << fn << endl;
//    // Read the file. If there is an error, report it and exit.

//    try
//    {
//        cfg.readFile(fn);
//    }
//    catch(const FileIOException &fioex)
//    {
//      std::cerr << "I/O error while reading file." << std::endl;
//      return(EXIT_FAILURE);
//    }
//    catch(const ParseException &pex)
//    {
//      std::cerr << "Parse error at " << pex.getFile() << ":" << pex.getLine()
//                << " - " << pex.getError() << std::endl;
//      return(EXIT_FAILURE);
//    }
//    int n;

//    try
//    {
//      n = cfg.lookup("n");

//    }
//    catch(const SettingNotFoundException &nfex)
//    {
//      cerr << "No 'name' setting in configuration file." << endl;
//    }
//    cout << n << " " << EXIT_SUCCESS << endl;
    return(0);

}
