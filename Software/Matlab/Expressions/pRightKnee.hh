/*
 * Automatically Generated from Mathematica.
 * Fri 26 Apr 2024 13:36:14 GMT-05:00
 */

#ifndef PRIGHTKNEE_HH
#define PRIGHTKNEE_HH

#ifdef MATLAB_MEX_FILE
// No need for external definitions
#else // MATLAB_MEX_FILE


#include "math2mat.hpp"
#include "mdefs.hpp"

namespace SymFunction
{

  void pRightKnee_raw(double *p_output1, const double *var1);

  inline void pRightKnee(Eigen::MatrixXd &p_output1, const Eigen::VectorXd &var1)
  {
    // Check
    // - Inputs
    assert_size_matrix(var1, 7, 1);

	
    // - Outputs
    assert_size_matrix(p_output1, 1, 3);


    // set zero the matrix
    p_output1.setZero();


    // Call Subroutine with raw data
    pRightKnee_raw(p_output1.data(), var1.data());
    }
  
  
}

#endif // MATLAB_MEX_FILE

#endif // PRIGHTKNEE_HH
