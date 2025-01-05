/*
 * Automatically Generated from Mathematica.
 * Fri 26 Apr 2024 13:36:08 GMT-05:00
 */

#ifndef DJ_LEGL_HH
#define DJ_LEGL_HH

#ifdef MATLAB_MEX_FILE
// No need for external definitions
#else // MATLAB_MEX_FILE


#include "math2mat.hpp"
#include "mdefs.hpp"

namespace SymFunction
{

  void dJ_LegL_raw(double *p_output1, const double *var1);

  inline void dJ_LegL(Eigen::MatrixXd &p_output1, const Eigen::VectorXd &var1)
  {
    // Check
    // - Inputs
    assert_size_matrix(var1, 14, 1);

	
    // - Outputs
    assert_size_matrix(p_output1, 2, 7);


    // set zero the matrix
    p_output1.setZero();


    // Call Subroutine with raw data
    dJ_LegL_raw(p_output1.data(), var1.data());
    }
  
  
}

#endif // MATLAB_MEX_FILE

#endif // DJ_LEGL_HH
