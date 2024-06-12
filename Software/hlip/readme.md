# hlip package
This package implements the hlip controller for the planar bipedal robot. 
## package structure
This package contains three layers of abstraction: 
1. HLIP class captures the essence of reduced-order HLIP model which plans the desired step length with reference step length
2. FiveLinkWalker_HLIP class is used to plan desired joint trajectories for the planar bipedal robot based on desired outputs
3. hardware_controller node takes in desired hardware inputs, instantiates a FiveLinkWalker_HLIP class to plan in 200Hz and outputs the desired joint positions and velocitites  