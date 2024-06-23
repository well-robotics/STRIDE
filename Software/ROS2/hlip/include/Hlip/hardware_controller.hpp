// #include <nav_msgs/msg/odometry.hpp>
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

// #include "Expression_hardware/COMPosition.hh"
// #include "Expression_hardware/COM_velocity.hh"
// #include "Expression_hardware/pLeftToe.hh"
// #include "Expression_hardware/vLeftToe.hh"
// #include "Expression_hardware/pRightToe.hh"
// #include "Expression_hardware/vRightToe.hh"
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
// #include "tf2/LinearMath/Matrix3x3.h"
// #include "tf2/LinearMath/Quaternion.h"
// #include <tf2_geometry_msgs/tf2_geometry_msgs.hpp>
#include "Eigen/Dense"
#include "FiveLinkWalker_HLIP.hpp"
#include "moving_average_filter.hpp"

using namespace rclcpp;
using namespace std::chrono_literals;
using namespace Eigen;
using namespace SymFunction;
using std::placeholders::_1;