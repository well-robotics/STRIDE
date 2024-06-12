/*
 * Automatically Generated from Mathematica.
 * Sun 9 Jun 2024 16:07:22 GMT-05:00
 */

#ifndef MMAT6_ROBOT_ASSEMBLY_V3_STRAIGHT_LEG_TPU_FOOT_HH
#define MMAT6_ROBOT_ASSEMBLY_V3_STRAIGHT_LEG_TPU_FOOT_HH

#ifdef MATLAB_MEX_FILE
// No need for external definitions
#else // MATLAB_MEX_FILE


#include "math2mat.hpp"
#include "mdefs.hpp"

namespace SymFunction
{

  void Mmat6_Robot_Assembly_v3_straight_leg_TPU_foot_raw(double *p_output1, const double *var1);

  inline void Mmat6_Robot_Assembly_v3_straight_leg_TPU_foot(Eigen::MatrixXd &p_output1, const Eigen::VectorXd &var1)
  {
    // Check
    // - Inputs
    assert_size_matrix(var1, 7, 1);

	
    // - Outputs
    assert_size_matrix(p_output1, 7, 7);


    // set zero the matrix
    p_output1.setZero();


    // Call Subroutine with raw data
    Mmat6_Robot_Assembly_v3_straight_leg_TPU_foot_raw(p_output1.data(), var1.data());
    }
  
  
}

#endif // MATLAB_MEX_FILE

#endif // MMAT6_ROBOT_ASSEMBLY_V3_STRAIGHT_LEG_TPU_FOOT_HH
