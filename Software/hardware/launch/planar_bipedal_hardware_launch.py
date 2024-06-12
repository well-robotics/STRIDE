from launch import LaunchDescription
from launch_ros.actions import Node 
from launch_ros.substitutions import FindPackageShare 
from ament_index_python.packages import get_package_share_path
from launch.actions import IncludeLaunchDescription
from launch.launch_description_sources import PythonLaunchDescriptionSource
from launch.substitutions import PathJoinSubstitution
import os
from ament_index_python.packages import get_package_share_directory

def generate_launch_description():
    control_yaml_file_name = "params/control_param.yaml"
    control_yaml_file = os.path.join(get_package_share_path("hlip"), control_yaml_file_name)
    return LaunchDescription([
        Node(
            package = 'hardware',
            #namespace = 'project',
            executable = 'imu_pub',
            name = 'imu_pub_node',
            prefix="nice -n -19 taskset -c 1",
            remappings=[
            ('imu_pub_node', 'remapped_imu_pub_node'),
            ]
        ),
        Node(
            package = 'hardware',
            #namespace = 'project',
            prefix="nice -n -19 taskset -c 1",
            executable = 'roll_velo_pub',
            name = 'roll_velo_pub_node'
        ),
        Node(
            package = 'hardware',
            #namespace = 'project',
            executable = 'contact_pub',
            prefix="nice -n -19 taskset -c 1",
            name = 'contact_pub_node',
        ),
        Node(
            package = 'hardware',
            #namespace = 'project',
            executable = 'joint_pub',
            prefix="nice -n -19 taskset -c 0",
            name = 'joint_inform_pub_node'
        ),
        #
        Node(
            package="hlip",
            name="hardware_controller",
            executable="hardware_controller",
            prefix="nice -n -19 taskset -c 3",
            output = "screen",
            # arguments=[("__log_level:=debug")],
            emulate_tty=True,
            parameters=[control_yaml_file]
        )
    ])