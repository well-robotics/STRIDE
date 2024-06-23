/*
 * Automatically Generated from Mathematica.
 * Sun 9 Jun 2024 16:09:55 GMT-05:00
 */

#ifndef PLEFTKNEE_HH
#define PLEFTKNEE_HH

#ifdef MATLAB_MEX_FILE
// No need for external definitions
#else // MATLAB_MEX_FILE


#include "math2mat.hpp"
#include "mdefs.hpp"

namespace SymFunction
{

  void pLeftKnee_raw(double *p_output1, const double *var1);

  inline void pLeftKnee(Eigen::MatrixXd &p_output1, const Eigen::VectorXd &var1)
  {
    // Check
    // - Inputs
    assert_size_matrix(var1, 7, 1);

	
    // - Outputs
    assert_size_matrix(p_output1, 1, 3);


    // set zero the matrix
    p_output1.setZero();


    // Call Subroutine with raw data
    pLeftKnee_raw(p_output1.data(), var1.data());
    }
  
  
}

#endif // MATLAB_MEX_FILE

#endif // PLEFTKNEE_HH
