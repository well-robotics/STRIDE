/*
 * @class: a ROS2 node for realizing H-LIP based walking
        It success the FiveLinkWalker class for planning the step 
        It sends desired output to lower level controller to follow the trajectory
 * @author: Yuhao Huang, Yicheng Zeng
*/

#include <rclcpp/rclcpp.hpp>
#include <rmw/types.h>

#include <sensor_msgs/msg/imu.hpp>
#include <sensor_msgs/msg/joint_state.hpp>
#include <std_msgs/msg/float64.hpp>
#include <std_msgs/msg/float32.hpp>

#include "communication/msg/contact_force.hpp"
#include "communication/msg/actuator_cmds.hpp"
#include "communication/msg/contact_state.hpp"
#include "communication/msg/motion_commands.hpp"
#include "communication/msg/joint_state.hpp"
#include "communication/msg/pitch_state.hpp"
#include "Expressions/COMPosition.hh"
#include "Expressions/COM_velocity.hh"
#include "Expressions/pLeftToe.hh"
#include "Expressions/vLeftToe.hh"
#include "Expressions/pRightToe.hh"
#include "Expressions/vRightToe.hh"

#include <chrono>
#include <cstdint>
#include <cstdio>
#include <string>
#include <cstring>
#include <iostream>
#include <memory>
#include <mutex>
#include <new>
#include <string>
#include <thread>
#include <type_traits>
#include <vector>
#include <random>
#include "rclcpp/rclcpp.hpp"
#include "Eigen/Dense"
#include "output_controller.hpp"
#include "moving_average_filter.hpp"
#include <deque>

using namespace rclcpp;
using namespace std::chrono_literals;
using namespace Eigen;
using namespace SymFunction;
using std::placeholders::_1;
