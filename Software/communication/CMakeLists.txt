cmake_minimum_required(VERSION 3.5)
project(communication)

# Default to C++14
 
set(CMAKE_CXX_STANDARD 14)


if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
add_compile_options(-Wall -Wextra -Wpedantic)
endif()





# if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
#   add_compile_options(-Wall -Wextra -Wpedantic)
# endif()

# find dependencies
find_package(ament_cmake REQUIRED)
find_package(std_msgs REQUIRED)
find_package(geometry_msgs REQUIRED)
find_package(sensor_msgs REQUIRED)
find_package(nav_msgs REQUIRED)

find_package(rosidl_default_generators REQUIRED)

find_package(rclcpp REQUIRED)

# find_package(ros_utilities REQUIRED)
# find_package(description REQUIRED)
# find_package(eom REQUIRED)
# find_package(control_common_toolbox REQUIRED)

include_directories(include)

# add_library(communication_node INTERFACE)


#  ament_target_dependencies(communication_node INTERFACE
#   rclcpp   
#   # ros_utilities
#   # description
#   # eom
#   nav_msgs
#   # control_common_toolbox
#   )

# target_include_directories(communication_node INTERFACE
#   "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
#   "$<INSTALL_INTERFACE:include>"
#   )
  
# ament_export_targets("export_communication_node" HAS_LIBRARY_TARGET)

# install(DIRECTORY include/${PROJECT_NAME}/
# DESTINATION include/${PROJECT_NAME}
# )

# install(TARGETS communication_node
# EXPORT "export_communication_node"
# ARCHIVE DESTINATION lib
# LIBRARY DESTINATION lib
# RUNTIME DESTINATION lib
# INCLUDES DESTINATION include
# )


ament_export_dependencies(
rosidl_default_runtime
rclcpp  
# eom
description
# ros_utilities
nav_msgs
# control_common_toolbox
)


file(GLOB_RECURSE source ${CMAKE_CURRENT_SOURCE_DIR}/msg/* ${CMAKE_CURRENT_SOURCE_DIR}/srv/* )
rosidl_generate_interfaces(${PROJECT_NAME}
        "msg/ActuatorCmds.msg"
        "msg/MotionCommands.msg"
        "msg/ContactState.msg"
        "msg/ContactForce.msg"
        "msg/JointState.msg"
        "msg/PitchState.msg"
        "srv/SimulationReset.srv"
        DEPENDENCIES ament_cmake std_msgs geometry_msgs sensor_msgs 
        )
##  this generate the target with using rosidl



ament_package()