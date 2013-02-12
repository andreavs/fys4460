#ifndef ATOM_H
#define ATOM_H

#include <string>
#include <armadillo>


class Atom
{
public:
    Atom();
    const arma::vec3& getVel(){return vel;}
    void setVel(const arma::vec3& newVel){vel = newVel;}
    const arma::vec3& getPos(){return pos;}
    void setPos(const arma::vec3& newPos){pos = newPos;}
    std::string getName(){return atomName;}
    void setName(std::string newName){atomName = newName;}
    void setSystemIndex(int sn){systemIndex = sn;}
    int getSystemIndex(){return systemIndex;}
    void setCellNumber(int cn){cellNumber = cn;}
    int getCellNumber(){return cellNumber;}




private:
    std::string atomName;
    arma::vec3 pos;
    arma::vec3 vel;
    int systemIndex;
    int cellNumber;


};
#endif // ATOM_H
