# Open-Source Planar Biped Robot STRIDE
![Hardware STRIDE](https://github.com/well-robotics/STRIDE/raw/main/Hardware_STRIDE.gif)
![Simulated STRIDE](https://github.com/well-robotics/STRIDE/raw/main/Sim_STRIDE.gif)
## Introduction
STRIDE is a low cost, easy-to-build yet versatile and robust planar bipedal robot platform for both educational and research purposes. The robot walks in a circular path constrained by a boom with a four-bar linkage. Except the basic robot and boom designs, a modular rough terrain setup and a quatitative measurable disturbance injection system are presented for system and algorithms validation purposes. The terrain setup can be easily reconfigured to fuilfill different research purposes while the disturbance injection system built with propellors provide measurable yet decoupled disturbances to test adaptation algorithms.  
This repo includes the software and hardware design of a open-sourced planar biped robot STRIDE and a step-to-step guid on how to make mechanical assemblies and how to implement a walking controller.
## Repository Structure
The package mainly consists of two folders, one for mechanical designs and one for software. The BOM includes all the materials needed to build the robot, along with an estimated cost for each item.
## Installation
### ROS2 
This project uses ROS2 as communication middleware. To configure ROS2 on workstation or PC, we recommend using Ubuntu 22.04 and ROS2 Humble which is validated by us. \
The installation guideline of ROS2 Humble is in: https://docs.ros.org/en/humble/Installation/Ubuntu-Install-Debians.html. To configure ROS2 on Raspberry Pi 4b+, you should build from source, follow the guidelines in https://docs.ros.org/en/humble/How-To-Guides/Installing-on-Raspberry-Pi.html and https://docs.ros.org/en/humble/Installation/Alternatives/Ubuntu-Development-Setup.html. \
After installation, you can verify the setup by checking environment variables: 
```bash
  printenv | grep -i ROS
```
Once ROS2 is setup successsfully, you should first build a ROS2 workspace
```bash
mkdir -p ~/ros2_ws/src
cd ~/ros2_ws
colcon build
```
Then clone the folder
```bash
git clone https://github.com/well-robotics/STRIDE.git
```
### LibSerial
This project uses LibSerial as the serial communication middleware, the detailed installation can be found in the Github page of LibSerial: https://github.com/crayzeewulf/libserial. Following the guideline of CMake based installations and use default installation folder. 

### FROST
FROST is used to generate dynamic libraries for the project. To install FROST on Windows, follow the guidelines in https://ayonga.github.io/frost-dev/pages/installation.html. Several things should be noticed: 
1. The version of Mathematica should be equal or lower than 12.3
2. Remember to change the Mathematica version in mathrun.m
### qpOASES
qpOASES is used to solve formulated QP problems when solving for desired joint torques. To install qpOASES on Windows, follow https://github.com/coin-or/qpOASES for guidelines. Several things should be noticed: 
1. The Matlab interface should be recompiled using Visual Studio 

## Launch Example
Once all the dependencies are well configured, first build the whole package: 
```bash
cd ~/ros2_ws
colcon build
```
Then, you can run the commands below to launch simulation and hardware nodes
### MuJoco Simulation
To run the MuJoco simulation, run
```bash
ros2 launch mujoco simulation_launch.py
```
Then a MuJoco simulation window will show up and the robot will be hang in the air initially. 
### Hardware Nodes
To run hardware nodes, run 
```bash
ros2 launch hardware hardware_launch.py
```
This launch file will successfully launch all the nodes when the Arduino Mega, two IMUs and two contact switches are connected to the Raspberry Pi. 
