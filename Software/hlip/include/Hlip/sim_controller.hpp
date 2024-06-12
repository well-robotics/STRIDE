#pragma once

#include <nav_msgs/msg/odometry.hpp>
#include <rclcpp/rclcpp.hpp>
#include <rmw/types.h>

#include <sensor_msgs/msg/imu.hpp>
#include <sensor_msgs/msg/joint_state.hpp>
#include <std_msgs/msg/float64.hpp>

// #include <tf2/LinearMath/Quaternion.h>
// #include <tf2_ros/transform_broadcaster.h>

#include "communication/msg/contact_force.hpp"
#include "communication/msg/actuator_cmds.hpp"
// #include "communication/srv/simulation_reset.hpp"
#include "communication/msg/contact_state.hpp"
#include "communication/msg/motion_commands.hpp"


#include "Expressions/COMPosition.hh"
#include "Expressions/COM_velocity.hh"
#include "Expressions/pLeftToe.hh"
#include "Expressions/vLeftToe.hh"
#include "Expressions/pRightToe.hh"
#include "Expressions/vRightToe.hh"

// // #include "Expressions/COMPosition.hh"
// // #include "Expressions/COM_velocity.hh"
// #include "Expressions/J_COMPosition.hh"
// // #include "Expressions/pLeftToe.hh"
// // #include "Expressions/vLeftToe.hh"
// #include "Expressions/J_leftToe.hh"
// // #include "Expressions/pRightToe.hh"
// // #include "Expressions/vRightToe.hh"
// #include "Expressions/J_rightToe.hh"
// #include "Expressions/pelvis_ori.hh"
// #include "Expressions/J_pelvis_ori.hh"


using namespace rclcpp;
using namespace std::chrono_literals;
using std::placeholders::_1;

