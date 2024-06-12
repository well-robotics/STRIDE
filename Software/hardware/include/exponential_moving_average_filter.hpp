#include "math.h"

class ExponentialMovingAverageFilter{
// This class implement the exponential moving average filter
// The idea of the filter is to use y(i) = a*x + (1-a) * y(i-1)
public: 
    void set_coeff(double coeff_input){
        this->coeff = coeff_input; 
    }
    double calc(double input){
        double output = coeff*input + (1-coeff)*output_rev;
        output_rev = output; 
        return output; 
    }

private: 
    double output_rev = 0;
    double coeff = 0; 
};