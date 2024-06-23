# hardware package
This package contains sensor reading and serial communication nodes, and the controller launch file.

## How to Use

To use this package, please use following 

To run any node, use the following command line:

```
# Check the nodes in this package (in your workspace/src folder)
source ../install/setup.bash 
ros2 pkg executables hardware

# run any node to chceck actuators or sensors:
ros2 run hardware contact_pub

# in a new command line, run this to see if the switches are functioning properlly
source ../install/setup.bash 
ros2 topic list
ros2 topic echo /hardware/ContactState
```

To launch the template controller, run the following line(make sure you have hardware setup and powered properly): 

```
ros2 launch hardware planar_bipedal_hardware_launch.py
```

## Package Structure

The following nodes are included in this package:

1. contact_sensor_node: for contact sensor reading
2. imu_node: for reading pelvis angle and angular velocity
3. roll_velocity_node: for reading x-direction velocity
4. serial_communication_joint_node: for serial communication with Arduino
