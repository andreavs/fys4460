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
    double* getPos(double ret[]);
    double* getVel(double ret[]);
    std::string moleculeType;
    double* pos;
    double* vel;
};

#endif // PARTICLE_H
