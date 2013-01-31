#ifndef PARTICLE_H
#define PARTICLE_H
#include <string>
#include <iostream>

class Particle
{
public:
    Particle(double x, double y, double z, double vx, double vy, double vz);
    void setMoleculeType(std::string molecule);
    void setPos(double x, double y, double z);
    void setVel(double vx, double vy, double vz);
    std::string getMoleculeType();
    void getPos(double*);
    void getVel(double*);
    std::string moleculeType;
private:
    double* pos;
    double* vel;
};

#endif // PARTICLE_H
