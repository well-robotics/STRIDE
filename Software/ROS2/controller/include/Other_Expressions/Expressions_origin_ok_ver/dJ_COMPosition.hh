/*
 * Automatically Generated from Mathematica.
 * Wed 20 Mar 2024 20:03:25 GMT-05:00
 */

#ifndef DJ_COMPOSITION_HH
#define DJ_COMPOSITION_HH

#ifdef MATLAB_MEX_FILE
// No need for external definitions
#else // MATLAB_MEX_FILE


#include "math2mat.hpp"
#include "mdefs.hpp"

namespace SymFunction
{

  void dJ_COMPosition_raw(double *p_output1, const double *var1,const double *var2);

  inline void dJ_COMPosition(Eigen::MatrixXd &p_output1, const Eigen::VectorXd &var1,const Eigen::VectorXd &var2)
  {
    // Check
    // - Inputs
    assert_size_matrix(var1, 7, 1);
    assert_size_matrix(var2, 7, 1);

	
    // - Outputs
    assert_size_matrix(p_output1, 3, 7);


    // set zero the matrix
    p_output1.setZero();


    // Call Subroutine with raw data
    dJ_COMPosition_raw(p_output1.data(), var1.data(),var2.data());
    }
  
  
}

#endif // MATLAB_MEX_FILE

#endif // DJ_COMPOSITION_HH
