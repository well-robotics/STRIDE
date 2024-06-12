#ifndef LIP_HPP
#define LIP_HPP
#include "Eigen/Dense" 

using namespace Eigen; 

class LIP{
public: 
    LIP(); 
    // parameters used in LIP model 


    // methods used in LIP model 
    Vector2d LIP_dynamics(double z, double t, double x0, double v0); 
    double LIP_orbital_energy(double z, double x, double v); 


};
#endif