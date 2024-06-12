/*
 * Automatically Generated from Mathematica.
 * Wed 20 Mar 2024 20:03:43 GMT-05:00
 */

#ifndef VLEGL_HH
#define VLEGL_HH

#ifdef MATLAB_MEX_FILE
// No need for external definitions
#else // MATLAB_MEX_FILE


#include "math2mat.hpp"
#include "mdefs.hpp"

namespace SymFunction
{

  void vLegL_raw(double *p_output1, const double *var1);

  inline void vLegL(Eigen::MatrixXd &p_output1, const Eigen::VectorXd &var1)
  {
    // Check
    // - Inputs
    assert_size_matrix(var1, 14, 1);

	
    // - Outputs
    assert_size_matrix(p_output1, 2, 1);


    // set zero the matrix
    p_output1.setZero();


    // Call Subroutine with raw data
    vLegL_raw(p_output1.data(), var1.data());
    }
  
  
}

#endif // MATLAB_MEX_FILE

#endif // VLEGL_HH
