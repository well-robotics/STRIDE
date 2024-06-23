#include <chrono>
#include <cstdint>
#include <cstdio>
#include <string>
#include <cstring>
#include <iostream>
#include <memory>
#include <mutex>
#include <new>
#include <string>
#include <thread>
#include <type_traits>
#include <vector>
#include "rclcpp/rclcpp.hpp"
#include "tf2/LinearMath/Matrix3x3.h"
#include "tf2/LinearMath/Quaternion.h"
#include <tf2_geometry_msgs/tf2_geometry_msgs.hpp>

#include "sim_controller.hpp"
#include "Eigen/Dense"
#include "FiveLinkWalker_HLIP.hpp"
// #include "mujoco_sim_helper.hpp"
// #include "MuJoCoMessageHandler.h"
// #include "ros_utilities/ros_utilities.hpp"

// #define MUJOCO_PLUGIN_DIR "mujoco_plugin"

using namespace std::chrono_literals;
using namespace Eigen;
using namespace SymFunction;

enum foot_placement_states{
    LEFT_FOOT = 1,
    RIGHT_FOOT = 2,
    BOTH_FEET = 3,
    FLOATING = 0
};
enum contact_debounce{
    NOT_CONTACT = 0,
    TEMP_CONTACT = 1,
    SOLID_CONTACT = 2,
};

namespace simulation_ns
{
    class HlipController: public rclcpp::Node
    {
    public:
        HlipController(): Node("hlip_controller") ,name_prefix("simulation/")
        {
            RCLCPP_INFO(this->get_logger(), "controller node start, building publishers, subscribers and callbacks...");
            /* This node initialize a h-lip controller,
            which includes a step generlizer
            also, reads sensor readings in mujoco simulator, and publish simulator motor control(position, vel)
            */

            //declare the actuator command structure, for publishin
            // struct ActuatorCmds
            // {
            //     double time = 0.0;
            //     std::vector<std::string> actuators_name;
            //     std::vector<float> kp;
            //     std::vector<float> pos_des;
            //     std::vector<float> kd;
            //     std::vector<float> vel_des;
            //     std::vector<float> torque_feedforward;
            // };

            //declare parameters, determined by the real robot model
            if(true){
            this->declare_parameter("orbit_type",1); // 1 means P1 and 2 means P2
            this->declare_parameter("vdes",0.1);
            this->declare_parameter("COM_height",0.25);
            this->declare_parameter("TSSP",0.5);
            this->declare_parameter("TDSP",0.);
            this->declare_parameter("step_size",0.05); // vdes * Tssp = step_size
            this->declare_parameter("swing_height",0.015);
            this->declare_parameter("torso_angle", 0.);
            this->declare_parameter("pelvis_init_height",0.485);

            this->declare_parameter("kp_hip",0.1);
            this->declare_parameter("kp_knee",0.1);
            this->declare_parameter("kd_hip",0.1);
            this->declare_parameter("kd_knee",0.1);


            this->declare_parameter("com_height_test",0.1);
            this->declare_parameter("x_swing_test",0.1);
            this->declare_parameter("z_swing_test",0.1);
            this->declare_parameter("pitch_test",0.1);
            this->com_height_test = this->get_parameter("com_height_test").as_double();
            this->x_swing_test = this->get_parameter("x_swing_test").as_double();
            this->z_swing_test = this->get_parameter("z_swing_test").as_double();
            this->pitch_test = this->get_parameter("pitch_test").as_double();

            // get all the parameters for the hlip controller
            this->orbit = this->get_parameter("orbit_type").as_int();
            this->vdes = this->get_parameter("vdes").as_double();
            this->COM_height = this->get_parameter("COM_height").as_double();
            this->TSSP = this->get_parameter("TSSP").as_double();
            this->TDSP = this->get_parameter("TDSP").as_double();
            this->step_size = this->get_parameter("step_size").as_double();
            this->swing_height = this->get_parameter("swing_height").as_double();
            this->torso_angle = this->get_parameter("torso_angle").as_double();
            this->orbit_init = (this->orbit==1) ? P1 : P2;

            this->pelvis_init_height = this->get_parameter("pelvis_init_height").as_double();
            this->kp_hip = this->get_parameter("kp_hip").as_double();
            this->kp_knee = this->get_parameter("kp_knee").as_double();
            this->kd_hip = this->get_parameter("kd_hip").as_double();
            this->kd_knee = this->get_parameter("kd_knee").as_double();
            }

            // quality of service
            auto qos = rclcpp::QoS(rclcpp::KeepLast(1), rmw_qos_profile_sensor_data);

            // publishers remap to /simulation/actuators_cmds
            actuators_cmds_publisher = this->create_publisher<communication::msg::ActuatorCmds>(
                name_prefix+"actuators_cmds", qos);
            motion_cmds_publisher = this->create_publisher<communication::msg::MotionCommands>(
                name_prefix+"motion_cmds", qos);
            actuators_cmds_timer = this->create_wall_timer(
                5ms, std::bind(&HlipController::actuators_cmds_callback,this));

            // // subscribers to joint encoder, imu and contact sensSor
            joint_state_subscription = this->create_subscription<sensor_msgs::msg::JointState>(
                name_prefix+"joint_states", qos, std::bind(&HlipController::joint_state_callback, this, std::placeholders::_1));
            imu_subscription = this->create_subscription<sensor_msgs::msg::Imu>(
                name_prefix+"imu_data",qos, std::bind( &HlipController::imu_callback, this, std::placeholders::_1 ));
            contact_subscription = this->create_subscription<communication::msg::ContactState>(
                name_prefix+"ContactState", qos, std::bind( &HlipController::contact_callback, this, std::placeholders::_1 ));
            odom_subscription = this->create_subscription<nav_msgs::msg::Odometry>(
                name_prefix+"odom", qos, std::bind( &HlipController::odom_callback, this, std::placeholders::_1 ));
            sim_time_subscription = this->create_subscription<std_msgs::msg::Float64>(
                name_prefix+"sim_time", qos, std::bind( &HlipController::sim_time_callback, this, std::placeholders::_1 ));

            //timer_callbacks
            timer_gen = this->create_wall_timer(10ms, std::bind(&HlipController::timer_callback_hlip, this));

            //This node do the following steps:
            //  1. Initialize a HLIP class, parameters are pre-defined in YAML and the model xml file.
            //  2. Read in data from the topics sent from MuJoCo
            //  3. Using hlip based controller to generate step size, then y_des(swing foot traj, COM traj)
            //  4. using IK to generate desired position and velocity of joints, then publish.

            //initiate some message data
            this->released = false;

            this->contact_ptr_->left_contact.data = false;
            this->contact_ptr_->right_contact.data = false;
            this->actuator_cmds.pos_des = {0.0, 0.0, 0.0, 0.0};
            this->actuator_cmds.vel_des = {0.0, 0.0, 0.0, 0.0};
            this->actuator_cmds.kp = {0.0, 0.0, 0.0, 0.0};
            this->actuator_cmds.kd = {0.0, 0.0, 0.0, 0.0};
            this->actuator_cmds.torque_feedforward = {0.0, 0.0, 0.0, 0.0};
            this->actuator_cmds.actuators_name = {"right_hip","right_knee","left_hip","left_knee"};

            this->walker_class.FiveLinkWalker_HLIP_initialize(
                this->orbit_init,this->vdes,this->COM_height,this->TSSP,this->TDSP,this->step_size,this->swing_height,this->torso_angle);

        }
    private:
    //attributes inherate from FiveLinkWalker
        FiveLinkWalker_HLIP walker_class;
        std::string name_prefix;
        int orbit = 1;//default orbit=P1
        double vdes=0.0;
        double COM_height=0.0;
        double TSSP=0.0;
        double TDSP=0.0;
        double step_size=0.0;
        double swing_height=0.0;
        double torso_angle=0.0;
        orbit_type orbit_init=P1;
        stance_leg contact_foot=left_leg;
        stance_leg previous_contact = right_leg; //1 for left_foot, 2 for right
        double contact_time=0.0;
        double time_stamp=0.000;
        double roll=0.0,pitch=0.0,yaw=0.0;
        double pelvis_init_height=0.0;
        foot_placement_states state_now=FLOATING;
        foot_placement_states state_prev=FLOATING;

        double time_now=0.0;
        double phase=0.0;
        bool released=false;
        bool start_receiving_message = false;
        double  com_height_test=0.0,x_swing_test=0.0,
                z_swing_test=0.0,pitch_test=0.0;

        int debounce_tick = 0;
        contact_debounce LeftContact = NOT_CONTACT;
        contact_debounce RightContact = NOT_CONTACT;

        float kp_hip = 5. , kp_knee = 3.0;
        float kd_hip = 1.5, kd_knee = 1.0;

        //the messages struct read
        communication::msg::ContactState contact;
        sensor_msgs::msg::JointState joint_states;
        sensor_msgs::msg::Imu IMU_data;
        //the message sent
        communication::msg::ActuatorCmds actuator_cmds;
        communication::msg::MotionCommands motion_cmds;
        nav_msgs::msg::Odometry odom;
        double sim_time=0.0;

        //ptrs to be used in subscription callback
        std::shared_ptr<sensor_msgs::msg::JointState> joint_state_ptr_  = std::make_shared<sensor_msgs::msg::JointState>();
        std::shared_ptr<sensor_msgs::msg::Imu>        imu_ptr_          = std::make_shared<sensor_msgs::msg::Imu>();
        std::shared_ptr<communication::msg::ContactState> contact_ptr_  = std::make_shared<communication::msg::ContactState>();
        std::shared_ptr<nav_msgs::msg::Odometry> odom_ptr_ = std::make_shared<nav_msgs::msg::Odometry>();


        //shared pointers for publishers, subscribers and timer callback.
        rclcpp::TimerBase::SharedPtr actuators_cmds_timer;
        rclcpp::TimerBase::SharedPtr timer_gen;
        rclcpp::Publisher<communication::msg::ActuatorCmds>::SharedPtr actuators_cmds_publisher;
        rclcpp::Publisher<communication::msg::MotionCommands>::SharedPtr motion_cmds_publisher;

        rclcpp::Subscription<sensor_msgs::msg::JointState>::SharedPtr joint_state_subscription;
        rclcpp::Subscription<sensor_msgs::msg::Imu>::SharedPtr imu_subscription;
        rclcpp::Subscription<communication::msg::ContactState>::SharedPtr contact_subscription;
        rclcpp::Subscription<nav_msgs::msg::Odometry>::SharedPtr odom_subscription;
        rclcpp::Subscription<std_msgs::msg::Float64>::SharedPtr sim_time_subscription;

        void actuators_cmds_callback(){ //publisher callback for actuators in form of ActuatorCmds msg
            // // this callback function publish the existing actuator command at a constant rate
            // actuators_cmds_publisher->publish(this->actuator_cmds);
            this->actuator_cmds.header.stamp = rclcpp::Clock().now();
            // // this->actuator_cmds.header.frame_id = "";
            // auto actuator_cmd_message= communication::msg::ActuatorCmds();
            this->actuator_cmds.kp={this->kp_hip,this->kp_knee,this->kp_hip,this->kp_knee};
            this->actuator_cmds.kd={this->kd_hip,this->kd_knee,this->kd_hip,this->kd_knee};
            this->motion_cmds.header.stamp = rclcpp::Clock().now();
            // auto motion_cmd_message = communication::msg::MotionCommands();

            this->actuators_cmds_publisher->publish(this->actuator_cmds);
            // actuators_cmds_publisher->publish(this->actuator_cmds);
            this->motion_cmds_publisher->publish(this->motion_cmds);
            // RCLCPP_INFO(this->get_logger(),"I am not sending anything");
        }

        void timer_callback_hlip(){
            // 0. Check timer functioning
            // Please note that controller time is using sim time, but the time stamps are aligned using RCLCPP time stamp
            RCLCPP_INFO(this->get_logger(),"simTime: %f",this->sim_time);
            if(this->sim_time<0.02){
                this->released = false;
            }

            // 1. Read in q - floating base and joints, emplace into joint_state_fl
            // If no message sent from ROS2, keep it in starting condition, and break the callback in case anything happens
            VectorXd joint_pos = VectorXd::Zero(7);
            VectorXd joint_vel = VectorXd::Zero(7);
            VectorXd joint_state_fl = VectorXd::Zero(14);
            if(this->joint_state_ptr_->position.empty()){
                // RCLCPP_WARN(this->get_logger(), "joint state a null pointer");
                joint_pos << 0,this->pelvis_init_height,0,0,0,0,0;
                joint_vel << 0,0,0,0,0,0,0;
                ////TODO: insert this configuration into actuator_cmds
                return;
            }else{
                // 0,1,2: floating base x, z, pitch; 3,4,5,6: right_hip, right_knee, left_hip, left_knee
                joint_pos(0) = this->odom_ptr_->pose.pose.position.x;
                joint_pos(1) = this->odom_ptr_->pose.pose.position.y+this->pelvis_init_height;//relative COM height to world origin
                // Because the hinge joint at pelvis is set, its reading can be just joint position, say joint_states[0]
                for(int i=2; i<7; i++){ joint_pos(i)=this->joint_state_ptr_->position[i-2];}
                // Reconfiguring the kinematics from frost
                joint_pos(3) = -joint_pos(3);
                joint_pos(4) = M_PI - joint_pos(4);
                joint_pos(5) = joint_pos(5);
                joint_pos(6) = joint_pos(6)- M_PI;
                //Velocity
                joint_vel(0) = this->odom_ptr_->twist.twist.linear.x;
                joint_vel(1) = this->odom_ptr_->twist.twist.linear.y;
                joint_vel(2) = this->imu_ptr_->angular_velocity.y;
                for(int i=3; i<7; i++){ joint_vel(i)= - this->joint_state_ptr_->velocity[i];} //THIS ONE IS NEGATIVED!!!!
                joint_vel(5) = - joint_vel(5);
                joint_vel(6) = - joint_vel(6);
            }
            joint_state_fl << joint_pos,joint_vel;

            // RCLCPP_INFO(this->get_logger(),"jnt pos: %f %f %f %f %f %f %f.",
            //     joint_pos(0),joint_pos(1),joint_pos(2),joint_pos(3),joint_pos(4),joint_pos(5),joint_pos(6));

            // 2. Release the robot, for details refer to mujoco/sim_helper
            if( this->sim_time>=0.5 &&  (this->contact_ptr_->left_contact.data || this->contact_ptr_->right_contact.data) && !this->released ){
                this->released = true;
                RCLCPP_WARN(this->get_logger(),"release the robot after %f seconds",this->sim_time);
            }

            //Determine stance leg.
            switch (this->contact_foot)
            {
            case left_leg:
                if( (this->walker_class.T >=  1.0 * this->walker_class.TSSP)  && (this->contact_ptr_->right_contact.data) ){//
                    this->previous_contact = left_leg;
                    this->contact_foot = right_leg;
                }
                break;
            case right_leg:
                if( (this->walker_class.T >=  1.0 * this->walker_class.TSSP) && (this->contact_ptr_->left_contact.data) ){//
                    this->previous_contact = right_leg;
                    this->contact_foot = left_leg;
                }
                break;
            default:
                break;
            }

            // 3.1 if the robot is not released, hold in air with 0 configuration
            if(!this->released){
                RCLCPP_WARN(this->get_logger(),"Not Released");
                Vector4d input_state_data = VectorXd::Zero(4); //(this->com_height_test,this->x_swing_test,this->z_swing_test,this->pitch_test);
                Vector4d d_input_state_data = VectorXd::Zero(4); //(0,0,0,0);
                
                // for fixed initial state
                this->walker_class.SSP_output_des <<this->com_height_test,this->x_swing_test,this->z_swing_test,this->pitch_test;
                this->walker_class.SSP_output_vdes << 0,0,0,0;

                // for sinesoidal wave test
                // this->phase = Pi/(0.5) * this->sim_time;
                // this->walker_class.SSP_output_des <<
                //         this->com_height_test,
                //         0.05 * cos( this->phase ),
                //         0.05 * (1-sin(this->phase)),
                //         0.0;
                // this->walker_class.SSP_output_vdes << 
                //         0,
                //         -0.05 * Pi/(0.5) * sin(this->phase),
                //         -0.05 * Pi/(0.5) * cos( this->phase ),
                //         0;

                MatrixXd COM_position_abs = MatrixXd::Zero(1,3);
                COMPosition(COM_position_abs, joint_pos);
                MatrixXd COMvelocity = MatrixXd::Zero(3,1);
                COM_velocity(COMvelocity,joint_pos,joint_vel);
                MatrixXd stance_leg_position = MatrixXd::Zero(2,1);
                MatrixXd swing_leg_position_abs = VectorXd::Zero(2,1);
                MatrixXd stance_leg_velocity = VectorXd::Zero(2,1);
                MatrixXd swing_leg_velocity = VectorXd::Zero(2,1);
                pRightToe(stance_leg_position,joint_pos);
                pLeftToe(swing_leg_position_abs,joint_pos);
                vRightToe(stance_leg_velocity, joint_state_fl);
                vLeftToe(swing_leg_velocity,joint_state_fl);
                input_state_data << (COM_position_abs(2)-stance_leg_position(1)),
                        (swing_leg_position_abs(0)-stance_leg_position(0)),
                        (swing_leg_position_abs(1)-stance_leg_position(1)),
                        joint_pos(2);
                d_input_state_data <<   COMvelocity(2) - stance_leg_velocity(1) ,
                        swing_leg_velocity(0) - stance_leg_velocity(0) ,
                        swing_leg_velocity(1) - stance_leg_velocity(1),
                        joint_vel(2);
                // Vector2d input_COM_state = Vector2d::Zero(2);
                // input_COM_state << (COM_position_abs(0) - stance_leg_position(0)), COMvelocity(0);
                // this->walker_class.FiveLinkWalker_HLIP_renew(joint_pos,joint_vel,
                // input_state_data,d_input_state_data,this->contact_foot,input_COM_state,this->sim_time);
                // RCLCPP_INFO(this->get_logger(), "try new IK");
                // this->walker_class.state_q << joint_pos(0), joint_pos(1), joint_pos(2), 0.6030, 2.3849, -.6119, -3.9071;
                this->walker_class.own_q << joint_pos.segment<5>(2);
                this->walker_class.FiveLinkWalker_HLIP_IK();

                // RCLCPP_INFO(this->get_logger(), "desired_motor_pos...w/%f %f %f %f",joint_pos_des[3],joint_pos_des[4],joint_pos_des[5],joint_pos_des[6]);
                // RCLCPP_INFO(this->get_logger(),"IK converged, hooray!");
                VectorXd joint_pos_des = this->walker_class.FiveLinkWalker_HLIP_IK_solver.getSolution_Position();
                VectorXd joint_vel_des = this->walker_class.FiveLinkWalker_HLIP_IK_solver.getSolution_Velocity();
                // for(long i=0; i<joint_pos_des.size();i++){
                //     wrap(joint_pos_des(i));
                //     wrap(joint_vel_des(i));
                //     // joint_pos_des(i) = wrap(joint_pos_des(i));
                //     // joint_vel_des(i) = wrap(joint_vel_des(i));
                // }
                this->motion_cmds.com_height_des = this->walker_class.SSP_output_des(0);
                this->motion_cmds.com_height_actual = COM_position_abs(2) - stance_leg_position(1);
                this->motion_cmds.x_des = this->walker_class.SSP_output_des(1);
                this->motion_cmds.x_actual = swing_leg_position_abs(0) - stance_leg_position(0);
                this->motion_cmds.z_swing_des = this->walker_class.SSP_output_des(2);
                this->motion_cmds.z_swing_actual = swing_leg_position_abs(1) - stance_leg_position(1);
                this->motion_cmds.pitch_des = 0;
                this->motion_cmds.pitch_vel_des = joint_pos(2);
                this->motion_cmds.pitch_vel_act = joint_vel(2);
                this->motion_cmds.vel_des.x =this->walker_class.SSP_output_vdes(1);
                this->motion_cmds.vel_des.z =this->walker_class.SSP_output_vdes(2);
                this->motion_cmds.vel_des.y =this->walker_class.SSP_output_vdes(0);//borrow y for COM_Z
                this->motion_cmds.contact_foot = this->contact_foot==left_leg?0:1;

                for(int i=0; i<4; i++){
                    this->actuator_cmds.pos_des[i] = joint_pos_des(i+1);
                    this->actuator_cmds.vel_des[i] = -joint_vel_des(i+1);  //THIS ONE IS NEGATIVED!!!!
                }
                this->actuator_cmds.vel_des[2] = -this->actuator_cmds.vel_des[2];
                this->actuator_cmds.vel_des[3] = -this->actuator_cmds.vel_des[3];
                this->actuator_cmds.pos_des[0] =    -    this->actuator_cmds.pos_des[0];
                this->actuator_cmds.pos_des[1] = M_PI - this->actuator_cmds.pos_des[1];
                this->actuator_cmds.pos_des[2] =    +  this->actuator_cmds.pos_des[2];
                this->actuator_cmds.pos_des[3] = M_PI + this->actuator_cmds.pos_des[3];

                // RCLCPP_INFO(this->get_logger(),"current COM_Height is : %f",
                //     COM_position_abs(2) - stance_leg_position(1));

                // this->actuator_cmds.vel_des = {0,0,0,0};
                // this->actuator_cmds.pos_des = {0,0,0,0};
                // // this->actuator_cmds.pos_des = {-0.6030, M_PI-2.3849, -.6119, 3.9071-M_PI};
                // this->actuator_cmds.pos_des = {-0.6030, 0.7566, -.6119, 0.7655};
                // this->actuator_cmds.pos_des = {-0.5030, 0.7566, -.5119, 0.7655};

                // // For PD testing
                // this->phase = Pi/(0.5) * this->sim_time;
                // RCLCPP_INFO(this->get_logger(),"outputing sine of: %f",fmod(this->phase,M_PI));
                // this->actuator_cmds.pos_des[0] = Pi/6 * sin( this->phase );// A sin( omega * t)
                // this->actuator_cmds.vel_des[0] = Pi/1 * Pi/6 * cos( this->phase );
                // this->actuator_cmds.pos_des[1] = 0;// A sin( omega * t)
                // this->actuator_cmds.vel_des[1] = 0;
                // this->actuator_cmds.pos_des[2] = -Pi/6 * sin( this->phase );// A sin( omega * t)
                // this->actuator_cmds.vel_des[2] = - Pi/1 * Pi/6 * cos( this->phase );
                // this->actuator_cmds.pos_des[3] = 0;
                // this->actuator_cmds.vel_des[3] = 0;//for PD testing

            }
            
            // 3.2 if the robot is released and touches the ground, start planning steps
            else if(this->released){
                MatrixXd stance_leg_position = MatrixXd::Zero(2,1);
                MatrixXd swing_leg_position_abs = VectorXd::Zero(2,1);
                MatrixXd stance_leg_velocity = VectorXd::Zero(2,1);
                MatrixXd swing_leg_velocity = VectorXd::Zero(2,1);
                if (this->contact_foot == right_leg){
                    pRightToe(stance_leg_position,joint_pos);
                    pLeftToe(swing_leg_position_abs,joint_pos);
                    vRightToe(stance_leg_velocity, joint_state_fl);
                    vLeftToe(swing_leg_velocity,joint_state_fl);
                }
                else if (this->contact_foot == left_leg){
                    pLeftToe(stance_leg_position,joint_pos);
                    pRightToe(swing_leg_position_abs,joint_pos);
                    vLeftToe(stance_leg_velocity, joint_state_fl);
                    vRightToe(swing_leg_velocity,joint_state_fl);
                }
                else{
                RCLCPP_WARN(this->get_logger(),"wrong stance leg info");
                }

                MatrixXd COM_position_abs = MatrixXd::Zero(1,3);
                COMPosition(COM_position_abs, joint_pos);
                MatrixXd COMvelocity = MatrixXd::Zero(3,1);
                COM_velocity(COMvelocity,joint_pos,joint_vel);
                COMvelocity(0)=joint_vel(0);
                // 2.1.3.3 put all the data into input states
                Vector4d input_state_data  = Vector4d::Zero(4);
                Vector4d d_input_state_data = Vector4d::Zero(4);
                Vector2d input_COM_state = Vector2d::Zero(2);
                input_state_data << (COM_position_abs(2)-stance_leg_position(1)),
                                    (swing_leg_position_abs(0)-stance_leg_position(0)),
                                    (swing_leg_position_abs(1)-stance_leg_position(1)),
                                    joint_pos(2);
                d_input_state_data <<   COMvelocity(2) - stance_leg_velocity(1) ,
                                        swing_leg_velocity(0) - stance_leg_velocity(0) ,
                                        swing_leg_velocity(1) - stance_leg_velocity(1),
                                        joint_vel(2);

                input_COM_state << (COM_position_abs(0) - stance_leg_position(0)), (COMvelocity(0) - stance_leg_velocity(0));

                // RCLCPP_INFO(this->get_logger(), "y:%.2f %.2f %.2f %.2f.", input_state_data(0), input_state_data(1), input_state_data(2), input_state_data(3));
                // RCLCPP_INFO(this->get_logger(), "com relative:%f", input_COM_state(0));
                // RCLCPP_INFO(this->get_logger(), "com relative:%.3f - %.3f ", COMvelocity(0),stance_leg_velocity(0));
                RCLCPP_INFO(this->get_logger(), "T:%f", this->walker_class.T);

                //put in all the values calculated into HLIP
                this->walker_class.FiveLinkWalker_HLIP_renew(joint_pos,joint_vel,
                input_state_data,d_input_state_data,this->contact_foot,input_COM_state,this->sim_time);

                // RCLCPP_INFO(this->get_logger(),"renewal successed");
                // RCLCPP_INFO(this->get_logger(),"using these to plan desired:  %f %f %f ",
                // (this->walker_class.SSP_output_start(0) - this->walker_class.COM_height),//real COM height
                // this->walker_class.step_real, // real step size by deadbeat control
                // this->walker_class.SSP_output_start(1));//staring swing x

                this->walker_class.FiveLinkWalker_HLIP_plan_desired();

                // RCLCPP_INFO(this->get_logger(),"planning successed");

                //try to publish the output of bezier curve calc.
                //TO be published: 1. z_com 2. swing_foot_x 3.swing_foot_z
                this->motion_cmds.com_height_des = this->walker_class.SSP_output_des(0);
                this->motion_cmds.com_height_actual = COM_position_abs(2) - stance_leg_position(1);
                this->motion_cmds.x_des = this->walker_class.SSP_output_des(1);
                this->motion_cmds.x_actual = swing_leg_position_abs(0) - stance_leg_position(0);
                this->motion_cmds.z_swing_des = this->walker_class.SSP_output_des(2);
                this->motion_cmds.z_swing_actual = swing_leg_position_abs(1) - stance_leg_position(1);
                this->motion_cmds.pitch_des = 0;
                this->motion_cmds.pitch_vel_des = joint_pos(2);
                this->motion_cmds.pitch_vel_act = joint_vel(2);
                this->motion_cmds.vel_des.x =this->walker_class.SSP_output_vdes(1);
                this->motion_cmds.vel_des.z =this->walker_class.SSP_output_vdes(2);
                this->motion_cmds.vel_des.y =this->walker_class.SSP_output_vdes(0);//borrow y for COM_Z
                this->motion_cmds.contact_foot = (this->contact_foot==left_leg)?0:1;

                this->walker_class.FiveLinkWalker_HLIP_IK();

                MatrixXd IK_Jacobian = MatrixXd::Zero(5,5);
                MatrixXd IK_y = MatrixXd::Zero(5,1);
                IK_Jacobian << this->walker_class.FiveLinkWalker_HLIP_J_output(this->walker_class.own_q);
                IK_y << this->walker_class.FiveLinkWalker_HLIP_y_output(this->walker_class.own_q);
                RCLCPP_INFO(this->get_logger(),"\n.%3f .%3f .%3f .%3f .%3f",
                IK_y(0,0),IK_y(1,0),IK_y(2,0),IK_y(3,0),IK_y(4,0)
                );
                // this->walker_class.Log2txt(IK_Jacobian,"test1");
                // this->walker_class.Log2txt(IK_y,"test1");

                if(this->walker_class.SSP_output_des.size() == 0){
                    RCLCPP_WARN(this->get_logger(), "Failed calculating IK output");
                }else{
                    // RCLCPP_INFO(this->get_logger(), "size_of SSP_output_vdes is %ld...",this->walker_class.SSP_output_des.size());
                }

                VectorXd joint_pos_des = this->walker_class.FiveLinkWalker_HLIP_IK_solver.getSolution_Position();
                VectorXd joint_vel_des = this->walker_class.FiveLinkWalker_HLIP_IK_solver.getSolution_Velocity();
                // for(long i=0; i<joint_pos_des.size();i++){
                //      wrap(joint_pos_des(i));
                //     wrap(joint_vel_des(i));
                // }

                for(int i=0; i<4; i++){
                    // this->actuator_cmds.pos_des[i] = this->walker_class.SSP_output_des[i] ;
                    // this->actuator_cmds.vel_des[i] = this->walker_class.SSP_output_vdes[i];
                    this->actuator_cmds.pos_des[i] = joint_pos_des(i+1);
                    this->actuator_cmds.vel_des[i] = -joint_vel_des(i+1);  //THIS ONE IS NEGATIVED!!!!
                }

                this->actuator_cmds.vel_des[2] = -this->actuator_cmds.vel_des[2];
                this->actuator_cmds.vel_des[3] = -this->actuator_cmds.vel_des[3];
                this->actuator_cmds.pos_des[0] =    -    this->actuator_cmds.pos_des[0];
                this->actuator_cmds.pos_des[1] = M_PI - this->actuator_cmds.pos_des[1];
                this->actuator_cmds.pos_des[2] =    +  this->actuator_cmds.pos_des[2];
                this->actuator_cmds.pos_des[3] = M_PI + this->actuator_cmds.pos_des[3];

                this->actuator_cmds.kp = {this->kp_hip,this->kp_knee,this->kp_hip,this->kp_knee};
                this->actuator_cmds.kd = {this->kd_hip,this->kd_knee,this->kd_hip,this->kd_knee};

                RCLCPP_INFO(this->get_logger(),"step planned");
            }

            // RCLCPP_INFO(this->get_logger(),"endTimerCallback");

        }

        void joint_state_callback(const sensor_msgs::msg::JointState::SharedPtr joint_state_msg){
            this->joint_states.header.frame_id = joint_state_msg->header.frame_id;
            this->joint_states.header.stamp = joint_state_msg->header.stamp;

            // this->joint_states.position.resize(joint_state_msg->name.size());
            this->joint_state_ptr_->name.resize(        joint_state_msg->name.size());
            this->joint_state_ptr_->position.resize(    joint_state_msg->position.size());
            this->joint_state_ptr_->velocity.resize(    joint_state_msg->velocity.size());
            for (size_t k =0; k<this->joint_state_ptr_->name.size(); k++){
                this->joint_state_ptr_->name[k]=joint_state_msg->name[k];
                this->joint_state_ptr_->position[k]=joint_state_msg->position[k];
                this->joint_state_ptr_->velocity[k]=joint_state_msg->velocity[k];
            }
            // for (size_t k =0; k<joint_state_msg->name.size(); k++){
            //     this->joint_states.position[k]=joint_state_msg->position[k];
            //     this->joint_states.velocity[k]=joint_state_msg->velocity[k];
            //     this->joint_states.name[k]=joint_state_msg->name[k];
            // }
            // RCLCPP_INFO(this->get_logger(),"%ld",this->joint_state_ptr_->position.size());
        }


        void imu_callback(const sensor_msgs::msg::Imu::SharedPtr imu_msg){
            ////TODO: read in IMU data
            // this->IMU_data = *imu_msg;
            this->imu_ptr_->header.frame_id = imu_msg->header.frame_id;
            this->imu_ptr_->header.stamp = imu_msg->header.stamp;

            this->imu_ptr_->angular_velocity.x = imu_msg->angular_velocity.x;
            this->imu_ptr_->angular_velocity.y = imu_msg->angular_velocity.y;
            this->imu_ptr_->angular_velocity.z = imu_msg->angular_velocity.z;

            this->imu_ptr_->linear_acceleration.x = imu_msg->linear_acceleration.x;
            this->imu_ptr_->linear_acceleration.y = imu_msg->linear_acceleration.y;
            this->imu_ptr_->linear_acceleration.z = imu_msg->linear_acceleration.z;

            this->imu_ptr_->orientation.w = imu_msg->orientation.w;
            this->imu_ptr_->orientation.x = imu_msg->orientation.x;
            this->imu_ptr_->orientation.y = imu_msg->orientation.y;
            this->imu_ptr_->orientation.z = imu_msg->orientation.z;
        }

        void contact_callback(const communication::msg::ContactState::SharedPtr contact_msg){
            ////TODO: read in the contact states
            // this->contact = *contact_msg;
            this->contact_ptr_->header.frame_id = contact_msg->header.frame_id;
            this->contact_ptr_->header.stamp = contact_msg->header.stamp;

            this->contact_ptr_->left_contact.data = contact_msg->left_contact.data;
            this->contact_ptr_->right_contact.data = contact_msg->right_contact.data;
        }

        void odom_callback(const nav_msgs::msg::Odometry::SharedPtr odom_msg){
            // this->odom = *odom_msg;
            this->odom_ptr_->header.frame_id = odom_msg->header.frame_id;
            this->odom_ptr_->header.stamp = odom_msg->header.stamp;
            this->odom_ptr_->pose.pose.position.x = odom_msg->pose.pose.position.x;
            this->odom_ptr_->pose.pose.position.y = odom_msg->pose.pose.position.y;
            this->odom_ptr_->pose.pose.position.z = odom_msg->pose.pose.position.z;
            this->odom_ptr_->pose.pose.orientation.w = odom_msg->pose.pose.orientation.w;
            this->odom_ptr_->pose.pose.orientation.x = odom_msg->pose.pose.orientation.x;
            this->odom_ptr_->pose.pose.orientation.y = odom_msg->pose.pose.orientation.y;
            this->odom_ptr_->pose.pose.orientation.z = odom_msg->pose.pose.orientation.z;
            this->odom_ptr_->twist.twist.linear.x = odom_msg->twist.twist.linear.x;
            this->odom_ptr_->twist.twist.linear.y = odom_msg->twist.twist.linear.y;
            this->odom_ptr_->twist.twist.linear.z = odom_msg->twist.twist.linear.z;
            this->odom_ptr_->twist.twist.angular.x = odom_msg->twist.twist.linear.x;
            this->odom_ptr_->twist.twist.angular.y = odom_msg->twist.twist.linear.y;
            this->odom_ptr_->twist.twist.angular.z = odom_msg->twist.twist.linear.z;
        }

        void sim_time_callback(const std_msgs::msg::Float64 sim_time_msg){
            this->sim_time = sim_time_msg.data;
        }

        // double wrap(double angle){
        void wrap(double &angle){
            // angle = fmod( angle , M_PI ) - M_PI_2;
            // RCLCPP_INFO(this->get_logger(),"%f ",angle);
            
            angle = fmod(angle+M_PI,2*M_PI)-M_PI;


            // angle = fmod(angle+M_PI,2*M_PI);
            // if(angle<0) {angle+=2*M_PI;}
            // // angle-=M_PI_2;
            // if(angle>=M_PI_2 && angle<3*M_PI_2){
            //     angle-=M_PI;
            // }else if(angle>=3*M_PI_2 && angle<2*M_PI){
            //     angle-=2*M_PI;
            // }
            // RCLCPP_INFO(this->get_logger(),"%f %f %f",ang1, ang2, angle);
            
            // return angle;
        }
    };
}

int main(int argc, const char **argv){
    rclcpp::init(argc,argv);

    //create a ros2 node object and make it spin continuously.
    rclcpp::Node::SharedPtr hlip_controller_node =
        std::make_shared<simulation_ns::HlipController>();
    cout<<"Create node object"<<endl;
    //spining control
    rclcpp::spin(hlip_controller_node);
    rclcpp::shutdown();
    return 0;
}