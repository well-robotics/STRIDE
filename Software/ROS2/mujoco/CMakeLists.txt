cmake_minimum_required(VERSION 3.5)
project(mujoco)

# find dependencies
find_package(ament_cmake REQUIRED)

list(APPEND ${PROJECT_NAME}_CONFIG_EXTRAS
        cmake/CheckAvxSupport.cmake
        cmake/FindOrFetch.cmake
        cmake/MujocoHarden.cmake
        cmake/MujocoLinkOptions.cmake
        cmake/MujocoMacOS.cmake
        cmake/SimulateDependencies.cmake
        cmake/SimulateOptions.cmake
        )

# find dependencies
find_package(rclcpp REQUIRED)
find_package(glfw3 REQUIRED)
find_package(sensor_msgs REQUIRED)
find_package(std_msgs REQUIRED)
find_package(geometry_msgs REQUIRED)
find_package(communication REQUIRED)
find_package(tf2 REQUIRED)
find_package(tf2_ros REQUIRED)
find_package(nav_msgs REQUIRED)
find_package(cv_bridge REQUIRED)

# BLASFEO Option
if(NOT TARGET mujoco)
  add_library(mujoco SHARED IMPORTED)
  set_target_properties(mujoco PROPERTIES IMPORTED_LOCATION ${CMAKE_CURRENT_SOURCE_DIR}/lib/libmujoco.so)
endif()

file(GLOB_RECURSE source ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cc)

add_executable(simulation ${source})

target_link_libraries(simulation mujoco glfw pthread rt ${CMAKE_DL_LIBS})

ament_target_dependencies(simulation 
sensor_msgs 
std_msgs
rclcpp 
communication
description
geometry_msgs 
tf2 tf2_ros nav_msgs cv_bridge
)

target_include_directories(simulation PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>)
        


install(TARGETS simulation
        DESTINATION lib/${PROJECT_NAME})

install(DIRECTORY lib 
DESTINATION lib/mujoco)
        ### install the launch folder to the destination
install(DIRECTORY launch
         DESTINATION share/${PROJECT_NAME})

ament_package()
