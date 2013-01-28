#include <iostream>
#include "particle.h"

using namespace std;
int plussen(int n)
{
 return n+1;
}



int main()
{
    cout << "Hello World!" << endl;
    Particle *hei = new Particle();
    int a = plussen(2);
    return 0;
}

