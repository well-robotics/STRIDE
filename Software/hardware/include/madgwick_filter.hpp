#include <math.h>

class Madgwick_filter{
public:
    //Madgwick_filter(double deltat, double gyromeserror, double b):dt(deltat),gyroMeasError(gyromeserror),beta(b){}; 
    void update(){
        float norm;  // vector norm
        float SEqDot_omega_1, SEqDot_omega_2, SEqDot_omega_3, SEqDot_omega_4; // quaternion derivative from gyroscopes elements
        float f_1, f_2, f_3;
        // objective function elements
        float J_11or24, J_12or23, J_13or22, J_14or21, J_32, J_33;
        float SEqHatDot_1, SEqHatDot_2, SEqHatDot_3, SEqHatDot_4;

        // Axulirary variables to avoid reapeated calcualtions
        float halfSEq_1 = 0.5f * SEq_1;
        float halfSEq_2 = 0.5f * SEq_2;
        float halfSEq_3 = 0.5f * SEq_3;
        float halfSEq_4 = 0.5f * SEq_4;
        float twoSEq_1 = 2.0f * SEq_1;
        float twoSEq_2 = 2.0f * SEq_2;
        float twoSEq_3 = 2.0f * SEq_3;

        // Normalise the accelerometer measurement
        norm = sqrt(a_x * a_x + a_y * a_y + a_z * a_z);
        a_x /= norm;
        a_y /= norm;
        a_z /= norm;
    
         // Compute the objective function and Jacobian
        f_1 = twoSEq_2 * SEq_4- twoSEq_1 * SEq_3- a_x;
        f_2 = twoSEq_1 * SEq_2 + twoSEq_3 * SEq_4- a_y;
        f_3 = 1.0f- twoSEq_2 * SEq_2- twoSEq_3 * SEq_3- a_z;
        J_11or24 = twoSEq_3;
        J_12or23 = 2.0f * SEq_4;
        J_13or22 = twoSEq_1;
        J_14or21 = twoSEq_2;
        J_32 = 2.0f * J_14or21;
        J_33 = 2.0f * J_11or24;
        
        // Compute the gradient (matrix multiplication)
        SEqHatDot_1 = J_14or21 * f_2- J_11or24 * f_1;
        SEqHatDot_2 = J_12or23 * f_1 + J_13or22 * f_2- J_32 * f_3;
        SEqHatDot_3 = J_12or23 * f_2- J_33 * f_3- J_13or22 * f_1;
        SEqHatDot_4 = J_14or21 * f_1 + J_11or24 * f_2;
        
        // Normalise the gradient
        norm = sqrt(SEqHatDot_1 * SEqHatDot_1 + SEqHatDot_2 * SEqHatDot_2 + SEqHatDot_3 * SEqHatDot_3 + SEqHatDot_4 * SEqHatDot_4);
        SEqHatDot_1 /= norm;
        SEqHatDot_2 /= norm;
        SEqHatDot_3 /= norm;
        SEqHatDot_4 /= norm;
        
        // Compute the quaternion derrivative measured by gyroscopes
        SEqDot_omega_1 =-halfSEq_2 * w_x- halfSEq_3 * w_y- halfSEq_4 * w_z;
        SEqDot_omega_2 = halfSEq_1 * w_x + halfSEq_3 * w_z- halfSEq_4 * w_y;
        SEqDot_omega_3 = halfSEq_1 * w_y- halfSEq_2 * w_z + halfSEq_4 * w_x;
        SEqDot_omega_4 = halfSEq_1 * w_z + halfSEq_2 * w_y- halfSEq_3 * w_x;
        
        // Compute then integrate the estimated quaternion derrivative
        SEq_1 += (SEqDot_omega_1- (beta * SEqHatDot_1)) * dt;
        SEq_2 += (SEqDot_omega_2- (beta * SEqHatDot_2)) * dt;
        SEq_3 += (SEqDot_omega_3- (beta * SEqHatDot_3)) * dt;
        SEq_4 += (SEqDot_omega_4- (beta * SEqHatDot_4)) * dt;

        // Normalise quaternion
        norm = sqrt(SEq_1 * SEq_1 + SEq_2 * SEq_2 + SEq_3 * SEq_3 + SEq_4 * SEq_4);
        SEq_1 /= norm;
        SEq_2 /= norm;
        SEq_3 /= norm;
        SEq_4 /= norm;
    }

    void quaternion_to_Euler(){
        this->pitch_angle= asin(-2 * SEq_2 * SEq_4 + 2 * SEq_1* SEq_3); // pitch
        this->roll_angle = atan2(2 * SEq_3 * SEq_4 + 2 * SEq_1 * SEq_2, -2 * SEq_2 * SEq_2 - 2 * SEq_3* SEq_3 + 1); // roll
        this->yaw_angle = atan2(-2 * SEq_2 * SEq_3 - 2 * SEq_1 * SEq_4, 2 * SEq_3 * SEq_3 + 2 * SEq_4 * SEq_4 - 1); 
    }

    void set_data(double acc_x,double acc_y,double acc_z,double gyro_x,double gyro_y,double gyro_z){
        this->a_x = acc_x;
        this->a_y = acc_y;
        this->a_z = acc_z;
        this->w_x = gyro_x;
        this->w_y = gyro_y;
        this->w_z = gyro_z; 
    }

    double get_pitch_angle(){
        return this->pitch_angle; 
    }
private: 
    double dt = 0.005; 
    double gyroMeasError = 1.5; 
    double beta = sqrt(3.0/4.0) * gyroMeasError; 
    const double PI = 3.1415926535; 

    float a_x, a_y, a_z; // accelerometer measurement 
    float w_x, w_y, w_z; // gyroscope measurement
    float SEq_1 = 1.0, SEq_2 = 0.0, SEq_3 = 0.0, SEq_4 = 0.0;  // orientation quaternion 
    float roll_angle, pitch_angle, yaw_angle; 

};