/*
 * @class: implementation of H-LIP based step planner
 * @author: Yuhao Huang, Yicheng Zeng
*/
#include "HLIP.hpp"
#include "Eigen/Dense" 
#include <iostream> 
#include "math.h" 

using namespace Eigen;

void HLIP::HLIP_S2S_dynamics(){
    // if (this-> HLIP_orbit == P1){
    //     HLIP_P1_step(); 
    // }
    // else if(this -> HLIP_orbit == P2){
    //     HLIP_P2_step(); 
    // }
    // this->HLIP_final_state = this->A_S2S*this->HLIP_final_state + this->B_S2S*this->step_size; 
    HLIP_A_matrix(); 
    HLIP_B_matrix();
    HLIP_deadbeat_gain();
    // MatrixXd A = this->A_S2S;
    // MatrixXd B = this->B_S2S;
    // this->Kdeadbeat(0) = (A(0,0)*A(1,0)*B(0) + A(1,0)*A(1,1)*B(0) - (A(0,0)*A(0,0))*B(1) - A(0,1)*A(1,0)*B(1))/
    //         (-A(1,0)*(B(0)*B(0)) + A(0,0)*B(0)*B(1) - A(1,1)*B(0)*B(1) + A(0,1)*(B(1)*B(1))); 
    // this->Kdeadbeat(1) = (A(0,1)*A(1,0)*B(0) +((A(1,1)*A(1,1))*B(0) - A(0,0)*A(0,1)*B(1) - A(0,1)*A(1,1)*B(1)))/
    //         (-A(1,0)*(B(0)*B(0)) + A(0,0)*B(0)*B(1) - A(1,1)*B(0)*B(1) + A(0,1)*(B(1)*B(1))); 
    // Eigen::EigenSolver<Eigen::MatrixXd> solver(this->A_S2S+this->B_S2S*this->Kdeadbeat.transpose());
    // Eigen::VectorXd eig_K = solver.eigenvalues().real();
    // std::cout<<"eigen of A+BK " << eig_K << std::endl;
    Vector2d C;
    C << 0,0;
    this->model_params.transpose() << this->A_S2S,this->B_S2S,C; 
    // std::cout <<"A " << this->A_S2S <<"\n B " << this->B_S2S << "model_params" <<this->model_params <<std::endl; 
}
void HLIP::HLIP_A_matrix(){
    // this->A_S2S(0,0) = cosh(this->TSSP*this->lambda);
    // this->A_S2S(0,1) = this->TDSP*cosh(this->TSSP*this->lambda)+1/this->lambda*sinh(this->TSSP*this->lambda); 
    // this->A_S2S(1,0) = this->lambda*sinh(this->TSSP*this->lambda); 
    // this->A_S2S(1,1) = cosh(this->TSSP*this->lambda)+this->lambda*this->TDSP*sinh(this->TSSP*this->lambda);
    MatrixXd A(2,2); 
    A << cosh(this->TSSP*this->lambda),1.0/this->lambda*sinh(this->TSSP*this->lambda)
        ,this->lambda*sinh(this->TSSP*this->lambda),cosh(this->TSSP*this->lambda); 
    MatrixXd temp(2,2); 
    temp << 1,this->TDSP,0,1;
    Eigen::EigenSolver<Eigen::MatrixXd> solver(A);
    Eigen::VectorXd eig_A = solver.eigenvalues().real();
    std::cout<< eig_A << std::endl;
    
    //Eigen::MatrixXd TSSP_matrix = Eigen::MatrixXd::Constant(A.rows(), A.cols(), TSSP);
    this->A_S2S = A*temp; 
}
void HLIP::HLIP_B_matrix(){
    this->B_S2S(0) = - cosh(this->TSSP*this->lambda);
    this->B_S2S(1) = - this->lambda*sinh(this->TSSP*this->lambda);
}


// update class parameters
void HLIP::HLIP_set_system_parameters(orbit_type orbit, double COM_height, double TSSP, double TDSP){
    this-> HLIP_orbit = orbit; 
    this-> COM_height = COM_height; 
    this-> TSSP = TSSP;
    this-> TDSP = TDSP; 
    this-> lambda = sqrt(9.81/COM_height); 
}
void HLIP::HLIP_update_desired_walking_parameters(double vdes, double step_size){
    this-> vdes = vdes; 
    this-> step_size =step_size; 
}
void HLIP::HLIP_initialization(orbit_type orbit, double COM_height, double TSSP, double TDSP, double vdes, double step_size){
    HLIP_set_system_parameters(orbit, COM_height,TSSP,TDSP); 
    HLIP_update_desired_walking_parameters(vdes, step_size); 
    // initialize the adaptive dynamics
    HLIP_S2S_dynamics(); 
}
void HLIP::HLIP_set_stance_leg(stance_leg leg_id){
    this->HLIP_stance_leg = leg_id; 
}


// plan desired step size
void HLIP::HLIP_orbital_slope(){
    if(this->HLIP_orbit == P1){
        this-> P1_orbital_slope = (this->lambda/tanh(this->TSSP/2*this->lambda)); //cos->cot
    }
    else if (this->HLIP_orbit == P2)
    {
        this-> P2_orbital_slope = (this->lambda*tanh(this->TSSP/2*this->lambda));
    }
    else{
        std::cout<<"wrong orbit type"<<std::endl;  
    }
}
void HLIP::HLIP_P1_step(){
    HLIP_orbital_slope(); 
    // double step_size_es = this->vdes*(this->TSSP+this->TDSP);  // step_size v_dT = u*
    this->HLIP_final_state(0) = step_size/(2+this->TDSP*this->P1_orbital_slope); 
    this->HLIP_final_state(1) = this->HLIP_final_state(0)*this->P1_orbital_slope; 
}
void HLIP::HLIP_P2_step(){
    // with left leg as the stance foot first
    HLIP_orbital_slope(); 
    double d2 = (pow(this->lambda,2)*pow((1/cosh(this->lambda/2*this->TSSP)),2)*(this->TSSP+this->TDSP)*this->vdes)/(pow(this->lambda,2)*this->TDSP+2*this->P2_orbital_slope); 
    Vector2d sol; 
    this->HLIP_final_state_left(0) = (-this->step_size - this->TDSP * d2) / (2 + this->TDSP*this->P2_orbital_slope);
    this->HLIP_final_state_left(1) = this->P2_orbital_slope * this->HLIP_final_state_left(0) + d2;
    double right_step_size = this->vdes*(this->TSSP+this->TDSP) - this->step_size; 
    // with right leg as the stance foot next 
    this->HLIP_final_state_right(0) = (right_step_size - this->TDSP * d2) / (2 + this->TDSP*this->P2_orbital_slope);
    this->HLIP_final_state_right(1) = this->P2_orbital_slope * this->HLIP_final_state_right(0) + d2;
}
void HLIP::HLIP_deadbeat_gain(){
    // this->Kdeadbeat <<1, this->TDSP+1/this->lambda*cosh(this->TSSP*this->lambda)/sinh(this->TSSP*this->lambda); 
    this->Kdeadbeat << 0.8,0.04;//0.5,0.064 ;
}
void HLIP::HLIP_P1_get_step_size(Vector2d state_real){
    HLIP_P1_step(); 
    // HLIP_deadbeat_gain(); 
    this->step_real = this->step_size + this->Kdeadbeat.dot((state_real - this->HLIP_final_state));
    // map step with respect to boom length to maintain speed linear velocity
    // if(this->HLIP_stance_leg == right_leg){
    //     // left foot swing
    //     this->step_real = (1378.0/1500)*this->step_real;
    // }
    // else if (this->HLIP_stance_leg == left_leg){
    //     // right foot swing
    //     this->step_real = (1622.0/1500)*this->step_real; 
    // }
    std::cerr << "step size:" << this->step_real << "=" 
    << this->step_size << "+["<<this->Kdeadbeat(0) << "," << this->Kdeadbeat(1)
    << "] * [(" << state_real(0) << "-" << this->HLIP_final_state(0) <<"), "
    << "(" << state_real(1) << "-" << this->HLIP_final_state(1)<<")]." << std::endl; 
}
void HLIP::HLIP_P2_get_step_size(Vector2d state_real){
    // be sure to feed in right step_des for right foot and left foot 
    HLIP_P2_step(); 
    HLIP_deadbeat_gain(); 
    if (this-> HLIP_stance_leg == left_leg){
        this->step_real = this->step_size + this->Kdeadbeat.dot((state_real - this->HLIP_final_state_left));  
    }
    else if(this->HLIP_stance_leg == right_leg){
        this->step_real = this->step_size + this->Kdeadbeat.dot((state_real - this->HLIP_final_state_right)); 
    }
}

/////////// Adaptive HLIP ///////////
void HLIP::HLIP_estimate_dynamics(Vector2d HLIP_state_prev, double u_prev, Vector2d HLIP_state_now){
    // estimate model params for the S2S dynamics of HLIP
    // input: previous HLIP_state, previous step size, HLIP_state_now, previous model params(model_params)
    // output: new model params

    // gain matrix is chosen as identity
    MatrixXd gain = 1.0*MatrixXd::Identity(4,4);
    VectorXd phi = VectorXd::Zero(4,1); 
    phi<<HLIP_state_prev(0), HLIP_state_prev(1),u_prev, 1;
    this->model_params = this->model_params +  gain*phi*((phi.transpose()*phi).inverse())*
        (((HLIP_state_now - this->model_params.transpose()*phi)).transpose());
}

void HLIP::HLIP_deadbeat_gain_adaptive(){ 
    // calculate deadbeat gain based on new model params
    // input: model params matrix
    // output: deadbeat gain

    // extract system dynamics from the model_params
    MatrixXd A = this->model_params.transpose().leftCols(2);
    //std::cout << "A " << A << std::endl;
    MatrixXd B = this->model_params.transpose().col(2);
    //std::cout << "model_params " << this->model_params.transpose()<<std::endl;
    //std::cout << "B " << B << std::endl;
    this->Kdeadbeat(0) = (A(0,0)*A(1,0)*B(0) + A(1,0)*A(1,1)*B(0) - (A(0,0)*A(0,0))*B(1) - A(0,1)*A(1,0)*B(1))/
            (-A(1,0)*(B(0)*B(0)) + A(0,0)*B(0)*B(1) - A(1,1)*B(0)*B(1) + A(0,1)*(B(1)*B(1))); 
    this->Kdeadbeat(1) = (A(0,1)*A(1,0)*B(0) +((A(1,1)*A(1,1))*B(0) - A(0,0)*A(0,1)*B(1) - A(0,1)*A(1,1)*B(1)))/
            (-A(1,0)*(B(0)*B(0)) + A(0,0)*B(0)*B(1) - A(1,1)*B(0)*B(1) + A(0,1)*(B(1)*B(1))); 
    //std::cout << "deabeat gain " << "[ " << this->Kdeadbeat(0) << " " << this->Kdeadbeat(1) << " ]" <<std::endl;
}