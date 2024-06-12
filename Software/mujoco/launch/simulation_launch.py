import os
from ament_index_python.packages import get_package_share_path
from launch import LaunchDescription
from launch.actions import DeclareLaunchArgument
from launch.substitutions import LaunchConfiguration
from launch_ros.actions import Node


def generate_launch_description():
    # xml_file_name = "model/xml/bucky.xml"
    # xml_file_name = "model/planarBiped/RobotAssembly.xml"
    # xml_file = os.path.join(get_package_share_path("description"), xml_file_name)
    # xml_file = os.path.join("/home/welllab/planar_biped_ws/src/Planar-Bipedal-Software-Repo/description/model/planarBiped/RobotAssembly.xml")
    
    # hardware_xml_file = os.path.join("/home/welllab/planar_biped_ws/src/Planar-Bipedal-Software-Repo/description/model/planarBiped/Robot_Assembly_straight_knee_feet.xml")
    # hardware_xml_file_name = "model/planarBiped/Robot_Assembly_straight_knee_feet.xml"
    hardware_xml_file_name = "model/planarBiped/Robot_TPU_foot.xml"
    hardware_xml_file = os.path.join(get_package_share_path("description"), hardware_xml_file_name)

    control_yaml_file_name = "params/control_param.yaml"
    control_yaml_file = os.path.join(get_package_share_path("hlip"), control_yaml_file_name)
    # control_yaml_file = os.path.join("/home/welllab/planar_biped_ws/src/Planar-Bipedal-Software-Repo/hlip/params/control_param.yaml")
    
    # estimation_yaml_file_name = "estimator_param.yaml"
    # estimation_yaml_file = os.path.join(get_package_share_path("estimation"), estimation_yaml_file_name)
    
    return LaunchDescription([
        Node(
            package="mujoco",
            executable="simulation",
            name="simulation_mujoco",
            output="screen",
            parameters=[{"model_file": hardware_xml_file}],
            emulate_tty=True,
            prefix="nice -n -19 taskset -c 1",
            # prefix=["gdbserver localhost:3000"],
            # arguments=[("__log_level:=debug")],            
            # parameters=[ {"model_file": xml_file}, control_yaml_file],
            ),
        # Node(
        #     package = "hlip",
        #     executable="sim_controller",
        #     name = "sim_controller", ## this name
        #     # prefix=["gdbserver localhost:3000"],
        #     # prefix="nice -n -19 taskset -c 4,5",
        #     output = "screen",
        #     # arguments=[("__log_level:=debug")],
        #     emulate_tty=True, # maybe useful for setting parameters here.
        #     parameters=[control_yaml_file] 
        # ),
        ### testing hardware controller
        Node(
            package = "hlip",
            executable="hardware_controller",
            name = "hardware_controller", ## this name
            # prefix=["gdbserver localhost:3000"],
            prefix="nice -n -19 taskset -c 4",
            output = "screen",
            # arguments=[("__log_level:=debug")],
            emulate_tty=True, # maybe useful for setting parameters here.
            parameters=[control_yaml_file] 
        ),
        Node(
            package = "hlip",
            executable="sim_hard_bridge",
            name = "sim_hard_bridge",
            prefix="nice -n -19 taskset -c 3",
            output = "screen",
            emulate_tty=True, # maybe useful for setting parameters here.
            parameters=[control_yaml_file] 
        ),
        ])
    """   
  export RMW_IMPLEMENTATION=rmw_cyclonedds_cpp
    """    