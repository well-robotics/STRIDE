/*
 * Automatically Generated from Mathematica.
 * Sun 9 Jun 2024 16:09:27 GMT-05:00
 */

#ifndef GE_VEC_ROBOT_ASSEMBLY_V3_STRAIGHT_LEG_TPU_FOOT_HH
#define GE_VEC_ROBOT_ASSEMBLY_V3_STRAIGHT_LEG_TPU_FOOT_HH

#ifdef MATLAB_MEX_FILE
// No need for external definitions
#else // MATLAB_MEX_FILE


#include "math2mat.hpp"
#include "mdefs.hpp"

namespace SymFunction
{

  void Ge_vec_Robot_Assembly_v3_straight_leg_TPU_foot_raw(double *p_output1, const double *var1,const double *var2);

  inline void Ge_vec_Robot_Assembly_v3_straight_leg_TPU_foot(Eigen::MatrixXd &p_output1, const Eigen::VectorXd &var1,const Eigen::VectorXd &var2)
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
    Ge_vec_Robot_Assembly_v3_straight_leg_TPU_foot_raw(p_output1.data(), var1.data(),var2.data());
    }
  
  
}

#endif // MATLAB_MEX_FILE

#endif // GE_VEC_ROBOT_ASSEMBLY_V3_STRAIGHT_LEG_TPU_FOOT_HH
