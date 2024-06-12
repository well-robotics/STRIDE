#pragma once 
#ifndef FIVELINKWALKER_HLIP
#define FIVELINKWALKER_HLIP

#include "HLIP.hpp"
#include "Eigen/Dense"
#include "Eigen/Sparse" 
#include "bezier_tools.hpp"
#include "IK_newton.hpp"
#include "IK_Gaussian.hpp"
#include "IK_LM.hpp"

using namespace Eigen; 
using namespace bezier_tools; 

class FiveLinkWalker_HLIP: public HLIP{
    // inherit from HLIP class
public:
    // parameters 
    VectorXd state_q = VectorXd::Zero(7);  // 7 state variables, should be updated frequently
    VectorXd state_dq = VectorXd::Zero(7);  // 7 velocity of state variables, should be updated frequently  
    VectorXd own_q = VectorXd::Zero(5);

    VectorXd initial_state = VectorXd::Zero(7); //  
    stance_leg HLIP_previous_stance_leg = right_leg;  // Left leg is the first stance leg which means the previous_stance_leg should be the right leg 

    // HLIP state
    Vector2d robot_state = Vector2d::Zero(2); 
    // the x position and velocity of the robot COM, used in HLIP to get desired step size, should be updated frequently
    // this one can be calculated using pCOMposition and vCOMposition
    // because they are with respect to the stance foot point, we should first calculate the position of the stance foot then calculate the COMposition and minus the stance foot position (set the first three state variables to 0)

    // walking outputs 
    // 1. z-com:calculated from pCOMposition    2. x-sw:calculated from pfoot    3. z-sw:calculated from pfoot   4.pitch-pelvis: read from IMU data  
    // 1. dz-com: calculated from vCOMpostion   2. dx-sw:calculated from vfoot   3. dz-sw:calculated from vfoot  4.dpitch-pelvis: read from IMU data
    VectorXd SSP_output_des = VectorXd::Zero(4);  // built using bezier-curve 
    VectorXd SSP_output_vdes = VectorXd::Zero(4); 
    VectorXd SSP_output_real = VectorXd::Zero(4);  // reading from the robot, should be updated frequently
    VectorXd SSP_output_vreal = VectorXd::Zero(4); 
    
    Vector2d COM_x_real;

    // parameter for bezier curves 
    VectorXd SSP_output_start = VectorXd::Zero(4); // state at the start of a SSP, for TDSP = 0, record each time change the stance foot, used for constructing Bezier curves
    double swing_height;  // the value
    double torso_angle; 
    double swing_height_neg = -0.005; 
    double T= 0;  // this is the time the robot has passed during SSP, renew to 0 if leg changes   
    double time_start = 0; 
    
    // parameters planned by bezier curves, should be fed into SSP_output_des
    double bezier_swing_foot_z=0; 
    double bezier_swing_foot_dz=0; 
    double bezier_swing_foot_x=0; 
    double bezier_swing_foot_dx=0; 
    double bezier_COM_z=0; 
    double bezier_COM_dz=0; 

    // parameter for outputing PD contorl
    // IK solver
    // 1. call the initialize method of the IK class 
    // 2. pass the function pointer of forward kinematics 
    // 3. pass the function pointer of inverse kinematics 
    // example: ik_leftleg_calib_solver.J_fun_pointer = [this](VectorXd q_leg){ return this->J_leftfoot_orientation(q_leg); };
    // IK_Newton_Solver FiveLinkWalker_HLIP_IK_solver; 
    IK_Gaussian_Solver FiveLinkWalker_HLIP_IK_solver;
    // IK_LM_Solver FiveLinkWalker_HLIP_IK_solver;

    // before walking
    // we set the initial pose for all states to be 0
    // then directly use HLIP to plan a trajectory
    void FiveLinkWalker_HLIP_initialize(orbit_type orbit, double vdes, double COM_height, double TSSP, double TDSP, double step_size, double swing_height,double torso_angle);

    
// update the states of the robots
    // 1. the stance leg, 
    // 2. the present time in SSP, if stance leg changes then renew
    // 3. the current walking outputs 
    // 4. the current COM state
    void FiveLinkWalker_HLIP_renew(VectorXd state_q_real, VectorXd state_dq_real, VectorXd output_state_real, VectorXd doutput_state_real, stance_leg stance_leg_input, Vector2d COM_x,double T);


// The control goals 
// 1. control the vertical center of mass 
// 2. swing foot periodically lift and strike the ground 
// 3. control the horizontal position of the foot to control step size
// 4. control the orientation of pelvis and swing foot to a constant 
// control flow
// 1. determine step size from HLIP
// 2. design the horizontal trajectory of swing foot 
// 3. design vertical COM position trajectory 
// 4. design vertical position of the swing foot
// 5. set orientation for swing foot and pelvis 

//////////////// trajectory building ////////////////
// 1. plan swing foot z trajectory 
    void FiveLinkWalker_HLIP_swing_foot_z(); 
// 2. plan swing foot x trajectory 
    void FiveLinkWalker_HLIP_swing_foot_x(); 
// 3. plan COM z trajectory
    void FiveLinkWalker_HLIP_COM_z(); 
//4. plan COM x trajectory (only for doing IK)
    void FiveLinkWalker_HLIP_COM_x();
// construct all the trajectory in one go
    void FiveLinkWalker_HLIP_plan_desired();


/////////////// joint level PD /////////////
// use IK to find desired joint angles and velocities
    void FiveLinkWalker_HLIP_IK(); 
    MatrixXd FiveLinkWalker_HLIP_J_output(VectorXd q);
    MatrixXd FiveLinkWalker_HLIP_y_output(VectorXd q);
// task-space PD for desired torque
    VectorXd FiveLinkWalker_HLIP_output_PD(VectorXd q_real, VectorXd dotq_real, VectorXd q_des, VectorXd dotq_des, MatrixXd Jacob, MatrixXd K_p, MatrixXd K_d);
    MatrixXd FiveLinkWalker_HLIP_construct_J(VectorXd q);

// gravity compnesation
    VectorXd FiveLinkWalker_HLIP_gravity_compensation(VectorXd q); 

// Added to log the jacobians
    void Log2txt(const MatrixXd matrix, std::string filename);

/////////////// QP based control ////////////
// TODO: find whether qpOASES is usable. 
};
#endif