/*
 * Automatically Generated from Mathematica.
 * Wed 20 Mar 2024 20:03:39 GMT-05:00
 */

#ifndef DJ_HIP_HH
#define DJ_HIP_HH

#ifdef MATLAB_MEX_FILE
// No need for external definitions
#else // MATLAB_MEX_FILE


#include "math2mat.hpp"
#include "mdefs.hpp"

namespace SymFunction
{

  void dJ_hip_raw(double *p_output1, const double *var1);

  inline void dJ_hip(Eigen::MatrixXd &p_output1, const Eigen::VectorXd &var1)
  {
    // Check
    // - Inputs
    assert_size_matrix(var1, 14, 1);

	
    // - Outputs
    assert_size_matrix(p_output1, 3, 7);


    // set zero the matrix
    p_output1.setZero();


    // Call Subroutine with raw data
    dJ_hip_raw(p_output1.data(), var1.data());
    }
  
  
}

#endif // MATLAB_MEX_FILE

#endif // DJ_HIP_HH
