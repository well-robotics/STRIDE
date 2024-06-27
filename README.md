# Open-sourced planar biped robot STRIDE
## Introduction
STRIDE is a low cost, easy-to-build yet versatile and robust planar bipedal robot platform for both educational and research purposes. The robot walks in a circular path constrained by a boom with a four-bar linkage. Except the basic robot and boom designs, a modular rough terrain setup and a quatitative measurable disturbance injection system are presented for system and algorithms validation purposes. The terrain setup can be easily reconfigured to fuilfill different research purposes while the disturbance injection system built with propellors provide measurable yet decoupled disturbances to test adaptation algorithms.  
This repo includes the software and hardware design of a open-sourced planar biped robot STRIDE and a step-to-step guid on how to make mechanical assemblies and how to implement a walking controller.
## repository structure
The package mainly consists of two folders, one for mechanical designs and one for software. The BOM includes all the materials needed to build the robot, along with an estimated cost for each item.
## Installation
### LibSerial
This project uses LibSerial as the serial communication middleware, the detailed installation can be found in the Github page of LibSerial: https://github.com/crayzeewulf/libserial. Following the guideline of CMake based installations. 
