# Open-sourced planar biped robot STRIDE
## Introduction
STRIDE is a low cost, easy-to-build yet versatile and robust planar bipedal robot platform for both educational and research purposes. The robot walks in a circular path constrained by a boom with a four-bar linkage. Except the basic robot and boom designs, a modular rough terrain setup and a quatitative measurable disturbance injection system are presented for system and algorithms validation purposes. The terrain setup can be easily reconfigured to fuilfill different research purposes while the disturbance injection system built with propellors provide measurable yet decoupled disturbances to test adaptation algorithms.  
This repo includes the software and hardware design of a open-sourced planar biped robot STRIDE and a step-to-step guid on how to make mechanical assemblies and how to implement a walking controller.
## repository structure
The package mainly consists of two folders, one for mechanical designs and one for software. The BOM includes all the materials needed to build the robot, along with an estimated cost for each item.
## Installation
### ROS2 
This project uses ROS2 as communication middleware. To configure ROS2 on workstation or PC, the validated setup is Ubuntu 22.04 and ROS2 Humble. The installation guideline of ROS2 Humble is in: https://docs.ros.org/en/humble/Installation/Ubuntu-Install-Debians.html. 
To configure ROS2 on Raspberry Pi 4b+, you should build from source, follow the guidelines in https://docs.ros.org/en/humble/How-To-Guides/Installing-on-Raspberry-Pi.html and https://docs.ros.org/en/humble/Installation/Alternatives/Ubuntu-Development-Setup.html. 
After installation, you can verify it by checking environment variables: 
```bash
  printenv | grep -i ROS
```
### LibSerial
This project uses LibSerial as the serial communication middleware, the detailed installation can be found in the Github page of LibSerial: https://github.com/crayzeewulf/libserial. Following the guideline of CMake based installations and use default installation folder. 
## Launch Example
Once all the dependencies are well configured, 
