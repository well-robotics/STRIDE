/*
 * Automatically Generated from Mathematica.
 * Wed 21 Feb 2024 16:41:31 GMT-06:00
 */

#ifndef MMAT1_ROBOT_ASSEMBLY_V3_HH
#define MMAT1_ROBOT_ASSEMBLY_V3_HH

#ifdef MATLAB_MEX_FILE
// No need for external definitions
#else // MATLAB_MEX_FILE


#include "math2mat.hpp"
#include "mdefs.hpp"

namespace SymFunction
{

  void Mmat1_Robot_Assembly_v3_raw(double *p_output1, const double *var1);

  inline void Mmat1_Robot_Assembly_v3(Eigen::MatrixXd &p_output1, const Eigen::VectorXd &var1)
  {
    // Check
    // - Inputs
    assert_size_matrix(var1, 7, 1);

	
    // - Outputs
    assert_size_matrix(p_output1, 7, 7);


    // set zero the matrix
    p_output1.setZero();


    // Call Subroutine with raw data
    Mmat1_Robot_Assembly_v3_raw(p_output1.data(), var1.data());
    }
  
  
}

#endif // MATLAB_MEX_FILE

#endif // MMAT1_ROBOT_ASSEMBLY_V3_HH
