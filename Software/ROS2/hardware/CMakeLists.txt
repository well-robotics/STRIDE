cmake_minimum_required(VERSION 3.8)

project(hardware)

if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  add_compile_options(-Wall -Wextra -Wpedantic)
endif()

# find dependencies
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
set(THREADS_HAVE_PTHREAD_ARG 1)

find_package(ament_cmake REQUIRED)
find_package(rclcpp REQUIRED)
find_package(sensor_msgs REQUIRED)
find_package(std_msgs REQUIRED)
find_package(communication REQUIRED)
find_package(Eigen3 3.4 REQUIRED NO_MODULE)
add_executable(joint_pub src/joint_comm_node.cc)
add_executable(imu_pub src/imu_node.cc)
add_executable(contact_pub src/contact_sensor_node.cc)
add_executable(roll_velo_pub src/roll_velocity_node.cc)

include_directories(include)
include_directories(src)



ament_target_dependencies(joint_pub rclcpp sensor_msgs communication)
ament_target_dependencies(imu_pub rclcpp sensor_msgs communication)
ament_target_dependencies(contact_pub rclcpp communication)
ament_target_dependencies(roll_velo_pub rclcpp std_msgs)


set(IMU_FILES include/imu_i2c.h
              include/imu_lib.hpp
              src/imu_i2c.c
              src/imu_lib.cpp
              include/madgwick_filter.hpp)          
add_library(IMU_lib ${IMU_FILES})
target_link_libraries(imu_pub IMU_lib)
target_link_libraries(roll_velo_pub IMU_lib)

find_library(WIRINGPI_LIBRARIES NAMES wiringPi)
target_link_libraries(contact_pub ${WIRINGPI_LIBRARIES})

install(TARGETS joint_pub DESTINATION lib/${PROJECT_NAME})
install(TARGETS imu_pub DESTINATION lib/${PROJECT_NAME})
install(TARGETS contact_pub DESTINATION lib/${PROJECT_NAME})
install(TARGETS roll_velo_pub DESTINATION lib/${PROJECT_NAME})
install(DIRECTORY 
        launch 
        DESTINATION share/${PROJECT_NAME}/
)
install(DIRECTORY include/
  DESTINATION include
)

ament_export_include_directories(include)


# for libserial
find_package(Threads REQUIRED)
find_package(PkgConfig)
pkg_check_modules(SERIAL libserial)


target_include_directories(joint_pub PRIVATE ${SERIAL_INCLUDE_DIRS})
target_link_libraries(joint_pub ${SERIAL_LDFLAGS} ${CMAKE_THREAD_LIBS_INIT})
set_property(TARGET joint_pub PROPERTY CXX_STANDARD 11)


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
