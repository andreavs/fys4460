#include <iostream>
#include "atom.h"
#include "cell.h"
#include <list>
using namespace std;

int main()
{
    list<Atom*> alist;
    Atom* atom1 = new Atom();
    alist.push_back(atom1);
    list<Atom*>::iterator i;
    i = alist.begin();
    cout << (*i)->getName() << endl;
    cout << "Hello World!" << endl;
    int j = 3;
    int k;
    k = 3+j%j;
    cout << k << endl;

    vector<Atom*> atomList;
    atomList.push_back(new Atom());
    cout << atomList[0]->getName() << endl;
    k = 3.15/4;
    cout << k << endl;
    return 0;
}

