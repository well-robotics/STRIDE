// custom code for using ros2

#pragma once

// #include <geometry_msgs/msg/transform_stamped.hpp>
#include <nav_msgs/msg/odometry.hpp>
#include <rclcpp/rclcpp.hpp>
#include <rmw/types.h>
// #include <sensor_msgs/msg/image.hpp>
#include <sensor_msgs/msg/imu.hpp>
#include <sensor_msgs/msg/joint_state.hpp>
#include <std_msgs/msg/float64.hpp>

#include <tf2/LinearMath/Quaternion.h>
#include <tf2_ros/transform_broadcaster.h>

#include "communication/msg/contact_force.hpp"
#include "communication/msg/actuator_cmds.hpp"
#include "communication/msg/contact_state.hpp"
#include "communication/srv/simulation_reset.hpp"
// #include "communication/msg/joint_state.hpp"
// #include "communication/msg/ContactForce.hpp"
// #include "communication/msg/ActuatorCmds.hpp"
// #include "communication/srv/SimulationReset.hpp"
#include <random>
#include "array_safety.h"
#include "simulate.h"

// #include "control_common_toolbox/data_logger.hpp"

using namespace rclcpp;

using namespace std::chrono_literals;

namespace mujoco_ros2
{
    namespace mj = ::mujoco;
    namespace mju = ::mujoco::sample_util;

    class MuJoCoMessageHandler : public rclcpp::Node
    {
    public:
        struct ActuatorCmds
        {
            double time = 0.0;
            std::vector<std::string> actuators_name;
            std::vector<float> kp;
            std::vector<float> pos_des;
            std::vector<float> kd;
            std::vector<float> vel_des;
            std::vector<float> torque_feedforward;
        };

        double simTime = 0;

        MuJoCoMessageHandler(mj::Simulate *sim);
        ~MuJoCoMessageHandler();

        std::shared_ptr<ActuatorCmds> get_actuator_cmds_ptr();

        double time_start = this->now().nanoseconds(); 
    private:
        // Data_Logger logger;

        // void init_logging()
        // {
        //     logger.init("mujocoSim");
        //     logger.add_data(&simTime, "mujoco_sim_Time");
        //     logger.add_data(&actuator_cmds_ptr_->time, "ros_time");
        // };

        void reset_callback(
            const std::shared_ptr<communication::srv::SimulationReset::Request> request,
            std::shared_ptr<communication::srv::SimulationReset::Response> response
            );

        void imu_callback();

        void odom_callback();

        void joint_callback();

        void sim_time_callback();

        void actuator_feedback_callback(); // only send the ones with encoders

        void contact_callback();

        void contact_state_callback();

        void actuator_cmd_callback(
            const communication::msg::ActuatorCmds::SharedPtr msg) const;

        void parameter_callback(const rclcpp::Parameter &);

        void drop_old_message();

        mj::Simulate *sim_;
        std::string name_prefix, model_param_name;
        std::vector<rclcpp::TimerBase::SharedPtr>
            timers_;

        rclcpp::Publisher<sensor_msgs::msg::Imu>::SharedPtr
            imu_publisher_;

        rclcpp::Publisher<sensor_msgs::msg::JointState>::SharedPtr
            joint_state_publisher_;

        rclcpp::Publisher<sensor_msgs::msg::JointState>::SharedPtr
            actuator_state_publisher_;

        rclcpp::Publisher<std_msgs::msg::Float64>::SharedPtr
            sim_time_publisher_;

        rclcpp::Publisher<communication::msg::ContactForce>::SharedPtr
            contact_sensor_publisher_;

        rclcpp::Publisher<nav_msgs::msg::Odometry>::SharedPtr
            odom_publisher_;

        //added Feb 06
        rclcpp::Publisher<communication::msg::ContactState>::SharedPtr
            contact_state_publisher_;

        //added April 02
        // rclcpp::Publisher<communication::msg::JointState>::SharedPtr
        //     ctrl_out_publisher_;

        rclcpp::Subscription<communication::msg::ActuatorCmds>::SharedPtr
            actuator_cmd_subscription_;

        rclcpp::Service<communication::srv::SimulationReset>::SharedPtr
            reset_service_;

        std::shared_ptr<rclcpp::ParameterEventHandler>
            param_subscriber_;

        std::shared_ptr<rclcpp::ParameterCallbackHandle>
            cb_handle_;

        std::shared_ptr<ActuatorCmds>
            actuator_cmds_ptr_;



        std::thread spin_thread; // is this ever been used?
    };

} // namespace mujoco_ros2
