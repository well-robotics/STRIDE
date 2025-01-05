/*
 * Automatically Generated from Mathematica.
 * Fri 26 Apr 2024 13:35:04 GMT-05:00
 */

#ifndef PELVIS_ROTVEL_HH
#define PELVIS_ROTVEL_HH

#ifdef MATLAB_MEX_FILE
// No need for external definitions
#else // MATLAB_MEX_FILE


#include "math2mat.hpp"
#include "mdefs.hpp"

namespace SymFunction
{

  void pelvis_rotvel_raw(double *p_output1, const double *var1,const double *var2);

  inline void pelvis_rotvel(Eigen::MatrixXd &p_output1, const Eigen::VectorXd &var1,const Eigen::VectorXd &var2)
  {
    // Check
    // - Inputs
    assert_size_matrix(var1, 7, 1);
    assert_size_matrix(var2, 7, 1);

	
    // - Outputs
    assert_size_matrix(p_output1, 1, 1);


    // set zero the matrix
    p_output1.setZero();


    // Call Subroutine with raw data
    pelvis_rotvel_raw(p_output1.data(), var1.data(),var2.data());
    }
  
  
}

#endif // MATLAB_MEX_FILE

#endif // PELVIS_ROTVEL_HH
