#include "LIP.hpp"
#include "Eigen/Dense"
#include "math.h"
using namespace Eigen; 

Vector2d LIP::LIP_dynamics(double z, double t, double x0, double v0){
    Vector2d sol; 
    double Tc = pow(z/9.81,0.5); 
    sol(0) =  x0*cosh(t/Tc)+Tc*v0*sinh(t/Tc);
    sol(1) =  x0/Tc*sinh(t/Tc)+v0*cosh(t/Tc);  
    return sol; 
}
double LIP::LIP_orbital_energy(double z, double x, double v){
    return (v/2-pow(x,2)/(2*z)*9.81); 
}