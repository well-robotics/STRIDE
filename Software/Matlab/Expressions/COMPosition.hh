/*
 * Automatically Generated from Mathematica.
 * Fri 26 Apr 2024 13:35:20 GMT-05:00
 */

#ifndef COMPOSITION_HH
#define COMPOSITION_HH

#ifdef MATLAB_MEX_FILE
// No need for external definitions
#else // MATLAB_MEX_FILE


#include "math2mat.hpp"
#include "mdefs.hpp"

namespace SymFunction
{

  void COMPosition_raw(double *p_output1, const double *var1);

  inline void COMPosition(Eigen::MatrixXd &p_output1, const Eigen::VectorXd &var1)
  {
    // Check
    // - Inputs
    assert_size_matrix(var1, 7, 1);

	
    // - Outputs
    assert_size_matrix(p_output1, 1, 3);


    // set zero the matrix
    p_output1.setZero();


    // Call Subroutine with raw data
    COMPosition_raw(p_output1.data(), var1.data());
    }
  
  
}

#endif // MATLAB_MEX_FILE

#endif // COMPOSITION_HH
