#include <string>
#include <iostream>
#include <iomanip>
#include <cstdlib>
//#include <libconfig.h++>

#include "configreader.h"

using namespace libconfig;
using namespace std;

ConfigReader::ConfigReader(){
    //config = new Config();
    cout << "politically coorrect" << endl;
}

ConfigReader::ConfigReader(string filename)
{
    try
    {
        cfg.readFile(filename.c_str());
    }
    catch(const FileIOException &fioex)
    {
      std::cerr << "I/O error while reading file." << std::endl;

    }
    catch(const ParseException &pex)
    {
      std::cerr << "Parse error at " << pex.getFile() << ":" << pex.getLine()
                << " - " << pex.getError() << std::endl;

    }

}
