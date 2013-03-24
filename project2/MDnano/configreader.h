#ifndef CONFIGREADER_H
#define CONFIGREADER_H

#include <libconfig.h++>
using namespace libconfig;

class ConfigReader
{
public:
    ConfigReader();
    ConfigReader(std::string filename);
    Config cfg;
};

#endif // CONFIGREADER_H
