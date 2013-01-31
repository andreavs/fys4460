#include "particle.h"

Particle::Particle(double x, double y, double z, double vx, double vy, double vz)
{
    pos = new double[3];
    vel = new double[3];
    pos[0] = x; pos[1] = y; pos[2] = z;
    vel[0] = vx; vel[1] = vy; vel[2] = vz;
    //std::cout << pos[0] << std::endl;


}



void Particle::setMoleculeType(std::string molecule){
    Particle::moleculeType = molecule;
}

void Particle::setPos(double x, double y, double z){
    Particle::pos[0] = x; Particle::pos[1] = y; Particle::pos[2] = z;
}

void Particle::setVel(double vx, double vy, double vz){
    Particle::vel[0] = vx; Particle::vel[1] = vy; Particle::vel[2] = vz;
}

std::string Particle::getMoleculeType(){
    return Particle::moleculeType;
}

void Particle::getPos(double* ret){
    for(int i=0;i<3;i++){
        ret[i] = pos[i];
    }
}

void Particle::getVel(double* ret){
    for(int i=0;i<3;i++){
        ret[i] = vel[i];
    }
}

