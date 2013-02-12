#include <iostream>
#include "system.h"
#include <string>
using namespace std;
using namespace arma;

int main()
{
    cout << "Hello World!" << endl;

    System *mySystem = new System();
    mySystem->runSimulation();
    return 0;
}

