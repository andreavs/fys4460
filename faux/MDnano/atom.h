#ifndef ATOM_H
#define ATOM_H

#include <string>
#include <armadillo>

class ConfigReader;

class Atom
{
public:
    Atom();
    Atom(ConfigReader *cfgReader);
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

    void setInitialPos(const arma::vec3& initPos){initialPos = initPos;}
    const arma::vec3& getInitialPos(){return initialPos;}
    void setRealPos(const arma::vec3& initPos){realPos = initPos;}
    const arma::vec3& getRealPos(){return realPos;}

    void setMass(double m){ mass = m;}
    double getMass(){return mass;}





private:
    std::string atomName;
    arma::vec3 pos;
    arma::vec3 vel;
    int systemIndex;
    int cellNumber;
    arma::vec3 initialPos;
    arma::vec3 realPos;
    double mass;


};
#endif // ATOM_H
