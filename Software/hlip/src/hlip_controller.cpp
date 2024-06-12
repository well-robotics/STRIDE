#include <chrono> 
#include <functional>
#include <memory>
#include <string>
#include <stdio.h>

#include "rclcpp/rclcpp.hpp"
#include "std_msgs/msg/string.hpp"

#include "Eigen/Dense"
#include "Expressions/COMPosition.hh"
#include "Hlip/HLIP.hpp"

// subscribe topic: /imu_data /joint_states
// publish topics: /