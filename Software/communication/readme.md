# Communication package
This package contains all the self-defined message and service interfaces being used in simulation and hardware.
## hardware interfaces  
1. ActuatorCmds.msg: for sending desired joint trajectories from Raspberry Pi to Arduino
2. ContactState.msg: for sending contact sensor states from contact sensor node to walking controller
3. JointStates.msg: for sending joint encoder readings from serial communication node to walking controller 
4. PitchState.msg: for sending pelvis angle and pelvis angular velocity to walking controller
4. MotionCommands.msg: for outputs states logging

## simulation interfaces
1. MotionCommands.msg: for outputs states logging
2. SimulationReset.srv: for MuJoco simulation reset
