# mujoco Package

This package is modified from the [open sourced ROS2 MuJoCo package](https://github.com/MindSpaceInc/Spot-MuJoCo-ROS2) from mindspace.inc.

It creates a simulation thread in MuJoCo with MJCF model file provided by description package. And established a node to send multiple sensor messages to ROS2.

## Package Structure

MujocoMessageHandler is the connection between MuJoCo and ROS2, it samples the required sensor data and publishes to according topics. This includes /ContactState, /JointState, /ImuData.(Also check communication/msg for details).

simulation executable(main.cc) establishes the physics thread in MuJoCo

### Things to notice

1. **The hold-release procedure**: It is realized with [sim_helper](src/Planar-Bipedal-Software-Repo/mujoco/include/mujoco_sim_helper.hpp). It first holds the floating base of the robot by giving a extremely large damping and stiffness, then slowly release it to the ground.
   1. The time variables are set locally.
   2. The fixed degrees of freedom are the ones on the floating base. They will remain a small value after it is released.
2. **Applying control to MuJoCo**: check apply_ctrl() function in main.cc. It realized a PD+feedforward controller of motors in MuJoCo. The data is directly obtained from shared pointer in Message Handler Node.
3. **Launch file**: The simulation-launch file launches simulator and controller simultaneously. The model xml file and the control yaml file are also specified here.