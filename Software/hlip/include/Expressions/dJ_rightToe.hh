/*
 * Automatically Generated from Mathematica.
 * Sun 9 Jun 2024 16:09:46 GMT-05:00
 */

#ifndef DJ_RIGHTTOE_HH
#define DJ_RIGHTTOE_HH

#ifdef MATLAB_MEX_FILE
// No need for external definitions
#else // MATLAB_MEX_FILE


#include "math2mat.hpp"
#include "mdefs.hpp"

namespace SymFunction
{

  void dJ_rightToe_raw(double *p_output1, const double *var1);

  inline void dJ_rightToe(Eigen::MatrixXd &p_output1, const Eigen::VectorXd &var1)
  {
    // Check
    // - Inputs
    assert_size_matrix(var1, 14, 1);

	
    // - Outputs
    assert_size_matrix(p_output1, 2, 7);


    // set zero the matrix
    p_output1.setZero();


    // Call Subroutine with raw data
    dJ_rightToe_raw(p_output1.data(), var1.data());
    }
  
  
}

#endif // MATLAB_MEX_FILE

#endif // DJ_RIGHTTOE_HH
