#include <functional>
#include <fstream>
#include "bezier_tools.hpp"
#include "Eigen/Dense"
#include "Eigen/Sparse" 
#include "unsupported/Eigen/MatrixFunctions"
#include "HLIP.hpp"
#include "FiveLinkWalker_HLIP.hpp"
#include "IK_newton.hpp" 
#include "Expressions/COMPosition.hh"
#include "Expressions/COM_velocity.hh"
#include "Expressions/J_COMPosition.hh"
#include "Expressions/pLeftToe.hh"
#include "Expressions/vLeftToe.hh"
#include "Expressions/J_leftToe.hh"
#include "Expressions/pRightToe.hh"
#include "Expressions/vRightToe.hh"
#include "Expressions/J_rightToe.hh"
#include "Expressions/pelvis_ori.hh"
#include "Expressions/J_pelvis_ori.hh"
#include "Expressions/Ge_vec_Robot_Assembly_v3_straight_leg.hh"
#include "Expressions/u_map_five_link_walker.hh"

using namespace Eigen;  
using namespace bezier_tools; 
using namespace SymFunction; 


void FiveLinkWalker_HLIP::FiveLinkWalker_HLIP_initialize(orbit_type orbit, double vdes, double COM_height, double TSSP, double TDSP, double step_size, double swing_height, double torso_angle){
    // update this function to satisfy the inheritence 
    HLIP_initialization(orbit, COM_height, TSSP, TDSP, vdes, step_size); 
    this->swing_height = swing_height;
    this->torso_angle = torso_angle; 
    this->HLIP_stance_leg = right_leg; 
    this->FiveLinkWalker_HLIP_IK_solver.initialize("DesJointPosition", 5, 5, 1e-8, 8, false);  // 1e-8 8
    this->FiveLinkWalker_HLIP_IK_solver.J_fun_pointer = [this](VectorXd q)
    {
        return this->FiveLinkWalker_HLIP_J_output(q); 
    }; 
    this->FiveLinkWalker_HLIP_IK_solver.y_fun_pointer = [this](VectorXd q)
    {
        return this->FiveLinkWalker_HLIP_y_output(q); 
    };  
}


void FiveLinkWalker_HLIP::FiveLinkWalker_HLIP_renew(VectorXd state_q_real, VectorXd state_dq_real, VectorXd output_state_real, VectorXd doutput_state_real, stance_leg stance_leg_input, Vector2d COM_state,double time_now){
    // this function renews the states of the robot required for generating control signals
    // 1. determine whether the stance leg changes
    if (stance_leg_input != this->HLIP_stance_leg ){
        this->HLIP_previous_stance_leg = this->HLIP_stance_leg; 
        this->HLIP_stance_leg = stance_leg_input; 
        this->SSP_output_start = output_state_real;
        this->time_start = time_now; 

        // record previous data for adaptive control
        Vector2d HLIP_state_prev = this->robot_state; 
        double step_prev = this->step_real;

        this->robot_state << COM_state;
        Vector2d xerr = Vector2d::Zero(2);
        xerr<<(COM_state(0) - output_state_real(1)), (COM_state(1)-doutput_state_real(1));

        // // Adaptive update system dynmaics and deadbeat gain 
        // std::cout << "for estimate dynamics: " << HLIP_state_prev.transpose() <<" "<< step_prev <<" "<< this->robot_state.transpose() << std::endl;
        if(step_prev!=0){
            HLIP_estimate_dynamics(HLIP_state_prev,step_prev,this->robot_state);
            HLIP_deadbeat_gain_adaptive(); 
        }


        // calculate step_des
        if(this->HLIP_orbit == P1){
            // HLIP_P1_get_step_size(this->robot_state);
            HLIP_P1_get_step_size(xerr);
        }
        else if(this->HLIP_orbit == P2){
            HLIP_P2_get_step_size(this->robot_state);
        }

        this->HLIP_start_state_real = this->robot_state;   
        this->HLIP_final_state_real = this->HLIP_final_state;
    }
    this->SSP_output_real = output_state_real;  // make the calculation from joint state to these values as a member function
    this->SSP_output_vreal = doutput_state_real; 
    this->robot_state = COM_state;
    this->state_q =  state_q_real; 
    this->own_q << this->state_q.segment<5>(2);
    // this->own_q(0)=0;
    // cout << "q in self coordinate: " << this->own_q.transpose() << endl;
    this->state_dq = state_dq_real; 

    // get the time passed during a single SSP
    this->T = time_now - this->time_start; // the time a SSP phase has passed between 0 - TSSP 
    if (this->T > this->TSSP){
        this->T = this->TSSP; 
    }

}

void FiveLinkWalker_HLIP::FiveLinkWalker_HLIP_swing_foot_z(){
    // this function calculates the desired swing foot z position at specific time t during ssp 
    // Input 
    // 1. start_swing_z position, 
    // 2. max_height
    // 3. the time T for this SPP phase 
    // Output:
    // the desired foot z position at that time as a double 
    VectorXd control_points = VectorXd::Zero(6); 
    control_points << this->SSP_output_start(2), this->swing_height / 3,
        this->swing_height, this->swing_height, this->swing_height / 3, // zsw_neg / 2,
        this->swing_height_neg;
    double t_bezier = this->T/this->TSSP; 
    this-> bezier_swing_foot_z = bezier(control_points,t_bezier); 
    this-> bezier_swing_foot_dz = dbezier(control_points,t_bezier); 
}
// overload swing_foot_z to get entire trajectory from lift off to touch down

void FiveLinkWalker_HLIP::FiveLinkWalker_HLIP_swing_foot_x(){
    // this function calculates the desired swing foot x position at specific time t during ssp
    // Input
    // 1. step length
    // 2. start state  
    //      start_state is the horizontal position of the swing foot w.r.t the stance foot at the beginning of the current SSP 
    // Output: 
    // the desired foot x position at that time as a double 
    VectorXd control_points = VectorXd::Zero(5); 
    control_points << 0,0,1,1,1; 
    double t_bezier = this->T/this->TSSP; 
    double bht = bezier(control_points, t_bezier); 
    double dbht = dbezier(control_points, t_bezier); 
     //step_real is the desired output of hlip, or the real step length
    this->bezier_swing_foot_x = ((1-bht)*this->SSP_output_start(1) + bht*this->step_real);

    // check the states used to plan x
    this->bezier_swing_foot_dx = -dbht*this->SSP_output_start(1) + dbht*this->step_real+ bht*this->vdes;                                                       
}

void FiveLinkWalker_HLIP::FiveLinkWalker_HLIP_COM_z(){
    // this function calculates the desired swing foot z position at specific time t during ssp
    // Input
    // 1. swing height
    // 2. start state  
    //      start_state is the vetical position of the swing foot w.r.t the stance foot at the beginning of the current SSP 
    // Output: 
    // the desired foot z position at that time as a double 
    VectorXd control_points = VectorXd::Zero(5); 
    control_points <<0,0,1,1,1; 
    double t_bezier = this->T/this->TSSP; 
    double bht = bezier(control_points, t_bezier); 
    double dbht = dbezier(control_points,t_bezier); 
    this-> bezier_COM_z = ((1-bht)*(this->SSP_output_start(0) - this->SSP_output_start(2) )+bht*this->COM_height); 
    this-> bezier_COM_dz = -dbht*(this->SSP_output_start(0) - this->SSP_output_start(2) )+dbht*this->COM_height;    
}

void FiveLinkWalker_HLIP::FiveLinkWalker_HLIP_COM_x(){
    // this function plans the desired COM-X trajectory
    // it is only used for IK

    // 1. plan the COM-X trajectory based on start state and desired end state with bezier curve
    // 1.1 calculate five knot points

    // calculate the estimated end state
    VectorXd knot1;
    VectorXd knot2; 
    VectorXd knot3;
    knot1 = 3/4*this->HLIP_start_state_real + 1/4*this->HLIP_final_state_real; 
    knot2 = 1/2*this->HLIP_start_state_real + 1/2*this->HLIP_final_state_real;
    knot3 = 1/4*this->HLIP_start_state_real + 3/4*this->HLIP_final_state_real;
    // 1.2 do bezier curve
    // std::cout << this->HLIP_final_state_real << std::endl;
    VectorXd control_points = VectorXd::Zero(5,1); 
    control_points << this->HLIP_start_state_real(0), knot1(0), knot2(0), knot3(0), this->HLIP_final_state_real(0); 
    double t_bezier = this->T/this->TSSP; 
    double bht = bezier(control_points, t_bezier); 
    double dbht = dbezier(control_points,t_bezier); 
    this->COM_x_real << bht, dbht;

}

void FiveLinkWalker_HLIP::FiveLinkWalker_HLIP_plan_desired(){
    FiveLinkWalker_HLIP_swing_foot_x(); 
    FiveLinkWalker_HLIP_swing_foot_z();
    FiveLinkWalker_HLIP_COM_z();
    FiveLinkWalker_HLIP_COM_x();
    this->SSP_output_des <<  this->bezier_COM_z,this->bezier_swing_foot_x,this->bezier_swing_foot_z,this->torso_angle; 
    this->SSP_output_vdes<<  this->bezier_COM_dz,this->bezier_swing_foot_dx,this->bezier_swing_foot_dz,0;
}   

void FiveLinkWalker_HLIP::FiveLinkWalker_HLIP_IK(){
    // call solve_from_current method of IK to get the desired joint q 
    // three input variables to pro's code 
    // 1. q is the initial guess
    // 2. Ydes is the desired output 
    // 3. dYdes is the desired output velocity 

    VectorXd y_des = VectorXd::Zero(5);
    VectorXd dy_des = VectorXd::Zero(5);

    ////////// use predefined values as input
    y_des << this->COM_x_real(0), this->SSP_output_des ;// 
    dy_des << this->COM_x_real(1), this->SSP_output_vdes;// 
    this->FiveLinkWalker_HLIP_IK_solver.solve_from_current(this->own_q, y_des, dy_des);
}

//NOTICE THAT BOTH J AND Y FUNCTI0N CALCULATES THE relative output
MatrixXd FiveLinkWalker_HLIP::FiveLinkWalker_HLIP_J_output(VectorXd q){ // q should be a 5x1 self_coord state
    // this function calculates the Jacobian matrix of output y
    // 1.x-com 2.z-com 3.x-sw 4. z-sw 5.pitch-pelvis  
    // combine the third row of J_COM,  J_Toe, J_pelvis_ori
    VectorXd q_all = VectorXd::Zero(7) ;
    q_all << 0,0,q.segment<5>(0);
    // q_all << 0,0,0,q.segment<4>(0);
    // q_all << 0,0,this->own_q;
    MatrixXd state_COM = MatrixXd::Zero(3,7); 
    J_COMPosition(state_COM, q_all); 
    MatrixXd pitch_pelvis = MatrixXd::Zero(1,7); 
    J_pelvis_ori(pitch_pelvis,q_all); 
    MatrixXd swing_toe = MatrixXd::Zero(2,7); 
    MatrixXd stance_toe = MatrixXd::Zero(2,7);
    if(this->HLIP_stance_leg == left_leg){ // right leg is the swing foot
        J_leftToe(stance_toe,q_all);
        J_rightToe(swing_toe,q_all); 
    }
    else if(this->HLIP_stance_leg == right_leg){
        J_leftToe(swing_toe,q_all); 
        J_rightToe(stance_toe,q_all);
    }
    MatrixXd J_output = MatrixXd::Zero(5,5); 
    J_output.block<1,5>(0,0) << state_COM.block<1,5>(0,2) - stance_toe.block<1,5>(0,2) ; //THIS IS WHERE THE ERROR COMES FROM.
    J_output.block<1,5>(1,0) << state_COM.block<1,5>(2,2) - stance_toe.block<1,5>(1,2) ;
    J_output.block<2,5>(2,0) << swing_toe.block<2,5>(0,2) - stance_toe.block<2,5>(0,2) ; 
    J_output.block<1,5>(4,0) << pitch_pelvis.block<1,5>(0,2); 
    return J_output; 
}

MatrixXd FiveLinkWalker_HLIP::FiveLinkWalker_HLIP_y_output(VectorXd q){
    // 1.x-com 2.z-com 3.x-sw 4. z-sw 5.pitch-pelvis  
    // combine the third row of COMPosition,  toe, pelvis_ori
    VectorXd q_all = VectorXd::Zero(7) ;
    q_all << 0.,0.,q.segment<5>(0);
    MatrixXd state_COM = MatrixXd::Zero(1, 3); 
    COMPosition(state_COM, q_all); 

    MatrixXd pitch_pelvis = MatrixXd::Zero(1, 1); 
    pelvis_ori(pitch_pelvis,q_all); 
    MatrixXd swing_toe = MatrixXd::Zero(2, 1); 
    MatrixXd stance_toe = MatrixXd::Zero(2, 1); 
    if(this->HLIP_stance_leg == left_leg){ // right leg is the swing foot
        pRightToe(swing_toe,q_all); 
        pLeftToe(stance_toe,q_all);
    }
    else if(this->HLIP_stance_leg == right_leg){
        pLeftToe(swing_toe,q_all); 
        pRightToe(stance_toe,q_all);
    }
    MatrixXd output = MatrixXd::Zero(5,1); 
    // cout << "need to know COMh:" << state_COM <<" stance toe xy"<< stance_toe.transpose() << endl;
    output.block<1,1>(0,0) = state_COM.block<1,1>(0, 0) - stance_toe.block<1,1>(0, 0); 
    output.block<1,1>(1,0) = state_COM.block<1,1>(0, 2) - stance_toe.block<1,1>(1, 0); 
    output.block<2,1>(2,0) = swing_toe.block<2,1>(0, 0) - stance_toe.block<2,1>(0, 0); 
    output.block<1,1>(4,0) = pitch_pelvis; 
    return output; 
}

void FiveLinkWalker_HLIP::Log2txt(const MatrixXd matrix, std::string filename)
{
    IOFormat CleanFmt(20, 0, ", ", "\n", "[", "]");
    std::string path = "/home//log_txt/";
    std::ofstream outfile(path + filename + ".txt");
    outfile << matrix.format(CleanFmt);
    outfile.close();
}

VectorXd FiveLinkWalker_HLIP::FiveLinkWalker_HLIP_output_PD(VectorXd q_real, VectorXd dotq_real, VectorXd q_des, VectorXd dotq_des, MatrixXd Jacob, MatrixXd K_p, MatrixXd K_d){
    return Jacob.transpose()*(K_p*(q_real - q_des)+K_d*(dotq_real - dotq_des)); 
}

MatrixXd FiveLinkWalker_HLIP::FiveLinkWalker_HLIP_construct_J(VectorXd q){
    // This function calculates the Jacobian needed for output PD
    // Input: 1*7 state vector
    // Output: 4*4 matrix
    VectorXd q_all = VectorXd::Zero(7) ;
    q_all << 0,0,q.segment<5>(0);
    MatrixXd state_COM = MatrixXd::Zero(3,7); 
    J_COMPosition(state_COM, q_all); 
    MatrixXd pitch_pelvis = MatrixXd::Zero(1,7); 
    J_pelvis_ori(pitch_pelvis,q_all); 
    MatrixXd swing_toe = MatrixXd::Zero(2,7); 
    MatrixXd stance_toe = MatrixXd::Zero(2,7);
    if(this->HLIP_stance_leg == left_leg){ // right leg is the swing foot
        J_leftToe(stance_toe,q_all);
        J_rightToe(swing_toe,q_all); 
    }
    else if(this->HLIP_stance_leg == right_leg){
        J_leftToe(swing_toe,q_all); 
        J_rightToe(stance_toe,q_all);
    }
    // TODO: carefully check the dimension and contents of all the Jacobian
    MatrixXd J_output = MatrixXd::Zero(4,4);
    J_output.block<1,4>(0,0) << state_COM.block<1,4>(2,3) - stance_toe.block<1,4>(1,3) ;
    J_output.block<2,4>(1,0) << swing_toe.block<2,4>(0,3) - stance_toe.block<2,4>(0,3) ; 
    J_output.block<1,4>(3,0) << 10.*pitch_pelvis.block<1,4>(0,3); 
    return J_output; 
}

VectorXd FiveLinkWalker_HLIP::FiveLinkWalker_HLIP_gravity_compensation(VectorXd q){
    // This methods computes torque outputs for compensating gravity
    // Input: 7x1 state 
    // output: 4x1 torque
    // TODO: 
    //    1. test whether sparse matrix is quicker or dense matrix is quicker
    //    2. set error tolerrance and iteration number to fulfil calculation criteria 

    // G(q) = J^T*F + Bu 
    // 1. construct J: jacobian, B: input choosing matrix, G(q) the gravity term
    MatrixXd contact_foot_jacob = MatrixXd::Zero(2,7);
    if(this->HLIP_stance_leg == left_leg){ 
        J_leftToe(contact_foot_jacob,q); 
    }
    else if(this->HLIP_stance_leg == right_leg){
        J_rightToe(contact_foot_jacob,q);
    }
    Eigen::SparseMatrix<double> contact_foot_jacob_sparse = contact_foot_jacob.sparseView();
    MatrixXd B = MatrixXd::Zero(7,4); 
    u_map_five_link_walker(B,q); // Gear ratio is needed for computing B, ours are 71.2
    Eigen::SparseMatrix<double> B_sparse = B.sparseView();
    
    MatrixXd G = VectorXd::Zero(7);
    Ge_vec_Robot_Assembly_v3_straight_leg(G,q,VectorXd::Zero(7));
    Eigen::SparseMatrix<double> G_sparse = G.sparseView();

    // 2. initialize the Eigen LeastSquareConjugateGradient solver 
    // 2.1 initialize matrix A
    SparseMatrix<double> A(7,6);
    A.leftCols(4) = B_sparse; 
    A.rightCols(2) = contact_foot_jacob_sparse.transpose(); 
    // 2.2 initialize solver
    LeastSquaresConjugateGradient<SparseMatrix<double>> solver; 
    solver.setTolerance(1e-12);
    solver.setMaxIterations(5);
    solver.compute(A); 
     
    // 3. use solver to solve the problem
    VectorXd output = solver.solve(G);
    return output;

}