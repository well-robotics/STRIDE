#ifndef EKF_HPP
#define EKF_HPP

// #include "rclcpp/rclcpp.hpp"
// #include <sensor_msgs/msg/imu.hpp>

#include <iostream>
#include <chrono>
#include <Eigen/Dense>
#include <Eigen/Geometry>

// tf broadcast
// #include "geometry_msgs/msg/transform_stamped.hpp"
// #include "geometry_msgs/msg/pose_stamped.hpp"

using namespace Eigen;

namespace ekf
{
    class EKF
    {
    private:
        /* data */
        // rclcpp::Subscription<sensor_msgs::msg::Imu>::SharedPtr imu_sub;
        // rclcpp::Subscription<geometry_msgs::msg::PoseStamped>::SharedPtr vo_pose_sub;
        // rclcpp::Publisher<sensor_msgs::msg::Imu>::SharedPtr publisher_filter_;

        // rclcpp::TimerBase::SharedPtr timer_;

        Vector3d accel_b_ = Vector3d::Zero();
        Vector3d angular_vel_b_ = Vector3d::Zero();

        int msg_num = 0;

        double dt_ = 0.005;//time interval

        bool init_vo_pose = false;
        bool init_filter = false;

        VectorXd vo_pose_quaternion_ = VectorXd::Zero(4);

        Vector3d gravity_ = Vector3d::Zero();
        VectorXd quternion_prior_;

        double time_0_;
        double imu_time_;
        double vo_time_;
        int discrete_time = 0;
        std::vector<Vector3d> angular_vel_stack_;
        std::vector<Vector3d> accel_stack_;
        std::vector<double> imu_time_stack_;
        std::vector<int> discrete_time_stack_;

        std::vector<Vector4d> filter_quaternion_stack_;
        std::vector<Matrix4d> Cov_stack_;

        bool vo_new_ = false;

    private:
        VectorXd quaternion_ = VectorXd::Zero(4);
        VectorXd quaternion_pred_ = VectorXd::Zero(4);
        VectorXd quaternion_correct_ = VectorXd::Zero(4);
        MatrixXd Cov_q_ = MatrixXd::Zero(4, 4);
        MatrixXd Cov_q_pred_ = MatrixXd::Zero(4, 4);
        MatrixXd Cov_q_correct_ = MatrixXd::Zero(4, 4);

        Matrix3d C_accel_ = Matrix3d::Zero();
        Matrix3d C_gyro_ = Matrix3d::Zero();
        Matrix4d C_vo_ = MatrixXd::Zero(4, 4);
    public:
        void EKF_init();
        // void raw_data_callback(const sensor_msgs::msg::Imu::SharedPtr msg);

        // void timerCallback();
        void update_raw_data(Vector3d linear_acceleration, Vector3d angular_velocity, double time_now);
        Quaternionf predict(const Eigen::Vector3d &angular_velocity);
        Quaternionf correct(const Vector3d &acceleration);
        void get_measurement();
        void gyro_nonlinear_predict(VectorXd &quaternion_pred,const VectorXd &quaternion,const Vector3d &gyro_reading, MatrixXd &Cov_q, MatrixXd &Cov_q_pred);
        void gyro_nonlinear_correct(VectorXd &quaternion_correct,const VectorXd &quaternion_pred,const Vector3d &accel_readings, MatrixXd &Cov_q_pred, MatrixXd &Cov_q_correct);
        
        void gyro_2_Ohm(const Vector3d &gyro_reading, MatrixXd &Ohm);
        void quat_2_W(const VectorXd &quaternion, MatrixXd &W);
        void quat_2_Rot(const VectorXd &quaternion, MatrixXd &Rotation);
        void quat_2_H(const VectorXd &quaternion, MatrixXd &H);
        void quat_norm(VectorXd &quaternion);
    };
 
}

#endif