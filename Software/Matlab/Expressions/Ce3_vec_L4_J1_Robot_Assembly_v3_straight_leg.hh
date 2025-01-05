/*
 * Automatically Generated from Mathematica.
 * Fri 26 Apr 2024 13:33:44 GMT-05:00
 */

#ifndef CE3_VEC_L4_J1_ROBOT_ASSEMBLY_V3_STRAIGHT_LEG_HH
#define CE3_VEC_L4_J1_ROBOT_ASSEMBLY_V3_STRAIGHT_LEG_HH

#ifdef MATLAB_MEX_FILE
// No need for external definitions
#else // MATLAB_MEX_FILE


#include "math2mat.hpp"
#include "mdefs.hpp"

namespace SymFunction
{

  void Ce3_vec_L4_J1_Robot_Assembly_v3_straight_leg_raw(double *p_output1, const double *var1,const double *var2);

  inline void Ce3_vec_L4_J1_Robot_Assembly_v3_straight_leg(Eigen::MatrixXd &p_output1, const Eigen::VectorXd &var1,const Eigen::VectorXd &var2)
  {
    // Check
    // - Inputs
    assert_size_matrix(var1, 7, 1);
    assert_size_matrix(var2, 7, 1);

	
    // - Outputs
    assert_size_matrix(p_output1, 7, 1);


    // set zero the matrix
    p_output1.setZero();


    // Call Subroutine with raw data
    Ce3_vec_L4_J1_Robot_Assembly_v3_straight_leg_raw(p_output1.data(), var1.data(),var2.data());
    }
  
  
}

#endif // MATLAB_MEX_FILE

#endif // CE3_VEC_L4_J1_ROBOT_ASSEMBLY_V3_STRAIGHT_LEG_HH
