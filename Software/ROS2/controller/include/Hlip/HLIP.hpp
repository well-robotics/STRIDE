/*
 * @class: implementation of H-LIP based step planner
 * @author: Yuhao Huang, Yicheng Zeng
*/
#ifndef HLIP_HPP
#define HLIP_HPP

#include "Eigen/Dense"

using namespace Eigen;
enum orbit_type {P1= 1, P2= 2};  // 1 for P1 orbit, 2 for P2 orbit
enum stance_leg {left_leg = 0, right_leg =1}; // left leg for 0, right leg for 1
class HLIP {
    /**************** utilities inside the class ***************************
    ***** parameters and variables *****/
    // 1. HLIP parameters: step_size, COM_height, forward_velocity 
    // 2. time parameters：SSP time: TSSP, DSP time: TDSP
    // 3. S2S dynamics parameters: A matrix and B matrix
    // 4. state variables x and control variables u
    // 5. controller parameters: deadbeat control, .... 
    /***** methods *****/
    // 1. initialize an object with given HLIP parameters (will time parameters be fixed?)
    // 2. a function calculate the S2S dynamics each time 
    // 3. a function renew the state variables each time 
    // 4. a function calculate the control variable each time
    // 5. a function calculate orbital energy
public: 
    // HLIP parameters 
    double TSSP; 
    double TDSP; 
    double step_size;  // with left leg as stance leg first
    double step_real = 0;  // the real step length for orbit stabilization
    double COM_height; 
    double vdes;
    double lambda; // sqrt(g/z0)  
    orbit_type HLIP_orbit; 
    stance_leg HLIP_stance_leg = right_leg; 
    // HLIP state 
    // P1 
    Vector2d HLIP_final_state = Vector2d::Zero(2);
    Vector2d HLIP_final_state_real = Vector2d::Zero(2); // the final state of HLIP with real step size after deadbeat control
    Vector2d HLIP_start_state_real = Vector2d::Zero(2);
    // P2
    Vector2d HLIP_final_state_left = Vector2d::Zero(2); 
    Vector2d HLIP_final_state_right = Vector2d::Zero(2); 
    // HLIP dynamics 
    MatrixXd A_S2S = MatrixXd::Zero(2,2);
    MatrixXd B_S2S = MatrixXd::Zero(2,1);  
    MatrixXd model_params = MatrixXd::Zero(4,2); // [A B C]^T 4*2 for adaptive control
    // orbital parameters 
    double P1_orbital_slope; 
    double P2_orbital_slope; 
    // control parameters 
    Vector2d Kdeadbeat; 
    /////////////////////////////// rebuild the functions to use class parameters rather than directly passing parameters //////////////////
    // methods 
    // update parameters and initialization 
    void HLIP_set_system_parameters(orbit_type orbit, double COM_height, double TSSP, double TDSP); 
    void HLIP_update_desired_walking_parameters(double vdes, double step_size); 
    void HLIP_initialization(orbit_type orbit, double COM_height, double TSSP, double TDSP, double vdes, double step_size); 
    void HLIP_set_stance_leg(stance_leg leg_id); 
    // S2S dynamics 
    void HLIP_A_matrix();
    void HLIP_B_matrix();
    void HLIP_S2S_dynamics(); 
    // Orbit calculations 
    void HLIP_orbital_slope();
    void HLIP_P1_step(); 
    void HLIP_P2_step();
    // Orbit stabilization 
    // through deadbeat control 
    void HLIP_deadbeat_gain();
    void HLIP_P1_get_step_size(Vector2d state_real); 
    void HLIP_P2_get_step_size(Vector2d state_real); 
    // Future: capture point control
    // Future: LQR control 

    // Adaptive control
    void HLIP_estimate_dynamics(Vector2d HLIP_state_prev, double u_prev, Vector2d HLIP_state_now); 
    void HLIP_deadbeat_gain_adaptive();
};
#endif