/*
 * Automatically Generated from Mathematica.
 * Wed 20 Mar 2024 20:03:06 GMT-05:00
 */

#ifndef PELVIS_ORI_HH
#define PELVIS_ORI_HH

#ifdef MATLAB_MEX_FILE
// No need for external definitions
#else // MATLAB_MEX_FILE


#include "math2mat.hpp"
#include "mdefs.hpp"

namespace SymFunction
{

  void pelvis_ori_raw(double *p_output1, const double *var1);

  inline void pelvis_ori(Eigen::MatrixXd &p_output1, const Eigen::VectorXd &var1)
  {
    // Check
    // - Inputs
    assert_size_matrix(var1, 7, 1);

	
    // - Outputs
    assert_size_matrix(p_output1, 1, 1);


    // set zero the matrix
    p_output1.setZero();


    // Call Subroutine with raw data
    pelvis_ori_raw(p_output1.data(), var1.data());
    }
  
  
}

#endif // MATLAB_MEX_FILE

#endif // PELVIS_ORI_HH
