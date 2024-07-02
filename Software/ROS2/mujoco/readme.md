# Mujoco Package

This package is modified from the [open sourced ROS2 MuJoCo package](https://github.com/MindSpaceInc/Spot-MuJoCo-ROS2) from mindspace.inc.

It creates a simulation thread in MuJoCo with MJCF model file provided by description package. And established a node to send multiple sensor messages to ROS2.

## Package Structure

`MujocoMessageHandler.cpp` is the connection between MuJoCo and ROS2, it samples the required sensor data and publishes to according topics. This includes `/ContactState`, `/JointState`, `/ImuData`.(Also check communication/msg for details).

simulation executable(main.cc) establishes the physics thread in MuJoCo

## Installation Guide

This guide referenced this [installation instructions][https://gist.github.com/saratrajput/60b1310fe9d9df664f9983b38b50d5da], To download MuJoCo, please go check its [releases][https://github.com/google-deepmind/mujoco/releases]. This version of STRIDE uses [2.3.2][https://github.com/google-deepmind/mujoco/releases/tag/2.3.2]. Take 2.3.2 as an example, after downloading the zip file:

```bash
cd /home/username/
mkdir .mujoco
tar -xvf mujoco232-linux-x86_64.tar.gz -C ~/.mujoco/
```

Add these lines to your `.bashrc`:

```bas
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/username/.mujoco/mujoco-2.3.2/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/lib/nvidia
export PATH="$LD_LIBRARY_PATH:$PATH"
export LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libGLEW.so
```

To test if the installation is successful, use the following line to test:

```bash
cd ~/.mujoco/mujoco210/bin
./simulate ../model/humanoid.xml
```

### Things to notice

1. **The hold-release procedure**: It is realized with [sim_helper](src/Planar-Bipedal-Software-Repo/mujoco/include/mujoco_sim_helper.hpp). It first holds the floating base of the robot by giving a extremely large damping and stiffness, then slowly release it to the ground.
   1. The time variables are set locally.
   2. The fixed degrees of freedom are the ones on the floating base. They will remain a small value after it is released.
2. **Applying control to MuJoCo**: check apply_ctrl() function in main.cc. It realized a PD+feedforward controller of motors in MuJoCo. The data is directly obtained from shared pointer in Message Handler Node.
3. **Launch file**: The simulation-launch file launches simulator and controller simultaneously. The model xml file and the control yaml file are also specified here.