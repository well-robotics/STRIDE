/*
 * Automatically Generated from Mathematica.
 * Fri 26 Apr 2024 13:35:17 GMT-05:00
 */

#ifndef DJ_PELVIS_POS_HH
#define DJ_PELVIS_POS_HH

#ifdef MATLAB_MEX_FILE
// No need for external definitions
#else // MATLAB_MEX_FILE


#include "math2mat.hpp"
#include "mdefs.hpp"

namespace SymFunction
{

  void dJ_pelvis_pos_raw(double *p_output1, const double *var1,const double *var2);

  inline void dJ_pelvis_pos(Eigen::MatrixXd &p_output1, const Eigen::VectorXd &var1,const Eigen::VectorXd &var2)
  {
    // Check
    // - Inputs
    assert_size_matrix(var1, 7, 1);
    assert_size_matrix(var2, 7, 1);

	
    // - Outputs
    assert_size_matrix(p_output1, 2, 7);


    // set zero the matrix
    p_output1.setZero();


    // Call Subroutine with raw data
    dJ_pelvis_pos_raw(p_output1.data(), var1.data(),var2.data());
    }
  
  
}

#endif // MATLAB_MEX_FILE

#endif // DJ_PELVIS_POS_HH
