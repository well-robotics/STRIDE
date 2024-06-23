import os
from ament_index_python.packages import get_package_share_path
from launch import LaunchDescription
from launch.actions import DeclareLaunchArgument
from launch.substitutions import LaunchConfiguration
from launch_ros.actions import Node


def generate_launch_description():
    hardware_xml_file_name = "model/planarBiped/STRIDE.xml"
    hardware_xml_file = os.path.join(get_package_share_path("description"), hardware_xml_file_name)

    control_yaml_file_name = "params/control_param.yaml"
    control_yaml_file = os.path.join(get_package_share_path("hlip"), control_yaml_file_name)
    
    return LaunchDescription([
        Node(
            package="mujoco",
            executable="simulation",
            name="simulation_mujoco",
            output="screen",
            parameters=[{"model_file": hardware_xml_file}],
            emulate_tty=True,
            prefix="nice -n -19 taskset -c 1",
            ),
        Node(
            package = "hlip",
            executable="hardware_controller",
            name = "hardware_controller",
            prefix="nice -n -19 taskset -c 4",
            output = "screen",
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