# Hlip package
This package implements the S2S dynamics based controller for STRIDE
## package structure
This package contains three layers of abstraction: 
1. HLIP class captures the essence of reduced-order HLIP model which plans the desired step size with reference step size
2. output_controller class is used to plan desired joint trajectories for the planar bipedal robot based on desired outputs
3. S2S_controller node takes in desired hardware inputs, instantiates a output_controller class to plan in 200Hz and outputs the desired joint positions and velocitites  