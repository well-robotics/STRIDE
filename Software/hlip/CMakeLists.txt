cmake_minimum_required(VERSION 3.8)
project(hlip)

if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  add_compile_options(-Wall -Wextra -Wpedantic)
endif()

# find dependencies
find_package(ament_cmake REQUIRED)
# uncomment the following section in order to fill in
# further dependencies manually.
# find_package(<dependency> REQUIRED)
find_package(rclcpp REQUIRED)
find_package(communication REQUIRED)
find_package(sensor_msgs REQUIRED)
find_package(std_msgs REQUIRED)
find_package(nav_msgs REQUIRED)
find_package(geometry_msgs REQUIRED)
find_package(Eigen3 3.4 REQUIRED NO_MODULE)
find_package(tf2 REQUIRED)
find_package(tf2_ros REQUIRED)
find_package(tf2_geometry_msgs REQUIRED)
find_package(hardware REQUIRED)

include_directories(
  include
  include/Hlip
  ${hardware_INCLUDE_DIRS}
  # ${sensor_msgs_LIBRARIES}
)
# message(SEND_ERROR "Variables: ${hardware_INCLUDE_DIRS}")
# include_directories(include/Hlip)

set(hlip_files include/Hlip/HLIP.hpp
               src/HLIP.cpp
               include/Hlip/LIP.hpp
               src/LIP.cpp
               include/Hlip/FiveLinkWalker_HLIP.hpp
               src/FiveLinkWalker_HLIP.cpp
               include/Hlip/bezier_tools.hpp
               src/bezier_tools.cpp
               include/Hlip/IK_newton.hpp
              #  include/Expression_hardware/COMPosition.cc
              #  include/Expression_hardware/COMPosition.hh
              #  include/Expression_hardware/pRightToe.cc
              #  include/Expression_hardware/pRightToe.hh
              #  include/Expression_hardware/pLeftToe.cc
              #  include/Expression_hardware/pLeftToe.hh
              #  include/Expression_hardware/J_COMPosition.cc
              #  include/Expression_hardware/J_COMPosition.hh
              #  include/Expression_hardware/J_leftToe.cc
              #  include/Expression_hardware/J_leftToe.hh
              #  include/Expression_hardware/J_rightToe.cc
              #  include/Expression_hardware/J_rightToe.hh
              #  include/Expression_hardware/pelvis_ori.cc
              #  include/Expression_hardware/J_pelvis_ori.cc
              #  include/Expression_hardware/pelvis_ori.hh
              #  include/Expression_hardware/J_pelvis_ori.hh
              #  include/Expression_hardware/vRightToe.hh
              #  include/Expression_hardware/vRightToe.cc
              #  include/Expression_hardware/vLeftToe.hh
              #  include/Expression_hardware/vLeftToe.cc
              #  include/Expression_hardware/COM_velocity.hh   
              #  include/Expression_hardware/COM_velocity.cc   
               include/Expressions/COMPosition.cc
               include/Expressions/COMPosition.hh
               include/Expressions/pRightToe.cc
               include/Expressions/pRightToe.hh
               include/Expressions/pLeftToe.cc
               include/Expressions/pLeftToe.hh
               include/Expressions/J_COMPosition.cc
               include/Expressions/J_COMPosition.hh
               include/Expressions/J_leftToe.cc
               include/Expressions/J_leftToe.hh
               include/Expressions/J_rightToe.cc
               include/Expressions/J_rightToe.hh
               include/Expressions/pelvis_ori.cc
               include/Expressions/J_pelvis_ori.cc
               include/Expressions/pelvis_ori.hh
               include/Expressions/J_pelvis_ori.hh
              #  include/Expressions/J_pelvis_pos.hh
              #  include/Expressions/J_pelvis_pos.cc
              #  include/Expressions/pelvis_pos.hh
              #  include/Expressions/pelvis_pos.cc
               include/Expressions/u_map_five_link_walker.hh
               include/Expressions/u_map_five_link_walker.cc
               include/Expressions/Ge_vec_Robot_Assembly_v3_straight_leg.hh
               include/Expressions/Ge_vec_Robot_Assembly_v3_straight_leg.cc
              )
add_library(Hlip_lib ${hlip_files})
target_link_libraries(Hlip_lib Eigen3::Eigen)
add_subdirectory(include/Expressions)
# add_subdirectory(include/Expression_hardware)

# add_executable(sim_controller src/simulator_controller.cpp)
# target_link_libraries(sim_controller Eigen3::Eigen)
# target_link_libraries(sim_controller DynamicExpressions)
# target_link_libraries(sim_controller Hlip_lib)

# ament_target_dependencies(sim_controller 
# sensor_msgs 
# nav_msgs
# std_msgs
# rclcpp 
# communication
# description
# geometry_msgs
# nav_msgs
# tf2_ros
# tf2
# tf2_geometry_msgs
# )

add_executable(hardware_controller src/hardware_controller.cpp)
target_link_libraries(hardware_controller Eigen3::Eigen)
target_link_libraries(hardware_controller DynamicExpressions)
target_link_libraries(hardware_controller Hlip_lib)

ament_target_dependencies(hardware_controller 
std_msgs
rclcpp 
communication
hardware
# description
)

# ament_target_dependencies(sim_controller rclcpp std_msgs) 
# install(TARGETS 
# sim_controller 
#   DESTINATION lib/${PROJECT_NAME})
ament_target_dependencies(hardware_controller rclcpp std_msgs hardware) 
install(TARGETS 
  hardware_controller 
  DESTINATION lib/${PROJECT_NAME})
install(FILES params/control_param.yaml
      DESTINATION share/${PROJECT_NAME}/params)

set(ekf_file src/ekf.cpp
             include/Hlip/ekf.hpp)
add_library(ekf_lib ${ekf_file})
target_link_libraries(ekf_lib Eigen3::Eigen)
# target_link_libraries(ekf_lib 
#   rclcpp
#   sensor_msgs 
#   ${sensor_msgs_LIBRARIES}
# )
# install(
#   TARGETS ekf_lib
#   ARCHIVE DESTINATION lib
#   LIBRARY DESTINATION lib
# )

# message(FATAL_ERROR "System processor: ${CMAKE_SYSTEM_PROCESSOR}")
# if (ARCHITECTURE_ID MATCHES "^(X86|X32|X64|x64|x86_64)$")
if (CMAKE_SYSTEM_PROCESSOR MATCHES "^(X86|X32|X64|x64|x86_64)$" )
  add_executable(sim_hard_bridge src/simulator_hardware_bridge.cpp)
  ament_target_dependencies(sim_hard_bridge 
  sensor_msgs 
  nav_msgs
  std_msgs
  rclcpp 
  communication
  # description
  geometry_msgs
  nav_msgs
  tf2_ros
  tf2
  tf2_geometry_msgs
  )
  target_link_libraries(sim_hard_bridge ekf_lib)
  install(
  TARGETS sim_hard_bridge 
  DESTINATION lib/${PROJECT_NAME})
endif()


if(BUILD_TESTING)
  find_package(ament_lint_auto REQUIRED)
  # the following line skips the linter which checks for copyrights
  # comment the line when a copyright and license is added to all source files
  set(ament_cmake_copyright_FOUND TRUE)
  # the following line skips cpplint (only works in a git repo)
  # comment the line when this package is in a git repo and when
  # a copyright and license is added to all source files
  set(ament_cmake_cpplint_FOUND TRUE)
  ament_lint_auto_find_test_dependencies()
endif()

ament_package()
