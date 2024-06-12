#include "rclcpp/rclcpp.hpp"
#include "hardware_controller.hpp"
#include <nav_msgs/msg/odometry.hpp>
#include "ekf.hpp"

using std::placeholders::_1;
using namespace std::chrono_literals;

class SimulatorToHardwarePublisher : public rclcpp::Node
{

  /* This class takes in simulator data and transfer them to what hardware controller can process to test hardware controller */

public:
  SimulatorToHardwarePublisher() : Node("simulator_to_hardware_publisher"), name_prefix("simulation/")
  {
    // setup the qos
    auto qos = rclcpp::QoS(rclcpp::KeepLast(1), rmw_qos_profile_sensor_data);
    //declare parameters
    declare_params();
    // initialize the publishers for publishing to topics listed:
    // 1. hardware/joint_states  communication::msg::JointState
    // 2. hardware/pitch         communication::msg::PitchState
    // 3. hardware/roll_velo     std_msgs::msg::Float32
    // 4. hardware/ContactState  communication::msg::ContactState
    pitch_publisher = this->create_publisher<communication::msg::PitchState>("hardware/pitch", qos);
    contact_publisher = this->create_publisher<communication::msg::ContactState>("hardware/ContactState", qos);
    joint_publisher = this->create_publisher<communication::msg::JointState>("hardware/joint_states", qos);
    roll_velo_publisher = this->create_publisher<std_msgs::msg::Float32>("hardware/roll_velo", qos);

    actuators_publihser = this->create_publisher<communication::msg::ActuatorCmds>("simulation/actuators_cmds",qos);
    motion_publisher = this->create_publisher<communication::msg::MotionCommands>("simulation/motion_cmds",qos);
    // actuators_cmds_timer = this->create_wall_timer(5ms, std::bind(&SimulatorToHardwarePublisher::commands_callback,this));

    // initialize timer for publishing
    timer_pub = this->create_wall_timer(5ms, std::bind(&SimulatorToHardwarePublisher::timer_callback, this));

    // initialize subscribers to the topics listed
    // 1.  simulator/joint_states  sensor_msgs::msg::JointState
    // 2.  simulator/imu_data      sensor_msgs::msg::Imu
    // 3.  simulator/ContactState  communication::msg::ContactState
    // 4.  simulator/odom          nav_msgs::msg::Odometry
    joint_state_subscription = this->create_subscription<sensor_msgs::msg::JointState>(
        name_prefix + "joint_states", qos, std::bind(&SimulatorToHardwarePublisher::joint_state_callback, this, std::placeholders::_1));
    imu_subscription = this->create_subscription<sensor_msgs::msg::Imu>(
        name_prefix + "imu_data", qos, std::bind(&SimulatorToHardwarePublisher::imu_callback, this, std::placeholders::_1));
    contact_subscription = this->create_subscription<communication::msg::ContactState>(
        name_prefix + "ContactState", qos, std::bind(&SimulatorToHardwarePublisher::contact_callback, this, std::placeholders::_1));
    odom_subscription = this->create_subscription<nav_msgs::msg::Odometry>(
        name_prefix + "odom", qos, std::bind(&SimulatorToHardwarePublisher::odom_callback, this, std::placeholders::_1));

    actuators_subscription = this->create_subscription<communication::msg::JointState>(
        "hardware/actuators_cmds", qos, std::bind(&SimulatorToHardwarePublisher::actuators_callback, this, std::placeholders::_1));
    motion_subscription = this->create_subscription<communication::msg::MotionCommands>(
        "hardware/motion_cmds", qos, std::bind(&SimulatorToHardwarePublisher::motion_callback, this, std::placeholders::_1));

    EKF_Filter.EKF_init();
  }

private:
  void timer_callback()
  {
    //1. extract positions and velocities from joint_states to joint_states_pub
    int i; 
    if(!this->joint_state_ptr_->position.empty()){
      for(i=0;i<4;i++){
        this->joint_states_pub.position[i] = this->joint_state_ptr_->position[i+1]; 
        this->joint_states_pub.velocity[i] = this->joint_state_ptr_->velocity[i+1]; 
      } 
    }
    joint_publisher->publish(this->joint_states_pub);

    //2. extract pitch angle and angular velocity from odom to pitch
    double imu_time = static_cast<double>(rclcpp::Clock().now().nanoseconds()) / 1e9 - this->time_start;
    // if(this->imu_ptr_->orientation.w!=0)
    // {
    //   Vector3d accel = Vector3d::Zero();
    //   Vector3d angular_vel = Vector3d::Zero();
    //   VectorXd quaternion = VectorXd::Zero(4);
    //   accel <<  this->imu_ptr_->linear_acceleration.x,
    //             this->imu_ptr_->linear_acceleration.y,
    //             this->imu_ptr_->linear_acceleration.z;
    //   angular_vel <<  this->imu_ptr_->angular_velocity.x, 
    //                   this->imu_ptr_->angular_velocity.y, 
    //                   this->imu_ptr_->angular_velocity.z;
    //   quaternion  <<  this->imu_ptr_->orientation.w,
    //                   this->imu_ptr_->orientation.x,
    //                   this->imu_ptr_->orientation.y,
    //                   this->imu_ptr_->orientation.z;
    //   
    //   // RCLCPP_INFO(this->get_logger(), "quat: %f, %f, %f, %f.",quaternion(0),quaternion(1),quaternion(2),quaternion(3));
    //   EKF_Filter.update_raw_data(accel,angular_vel,imu_time);
    //   Quaternionf q = EKF_Filter.predict(angular_vel);
    //   // RCLCPP_INFO(this->get_logger(), "quat: %f, %f, %f, %f.",quaternion(0),quaternion(1),quaternion(2),quaternion(3));
    //   q = EKF_Filter.correct(accel);

    //   auto euler = q.toRotationMatrix().eulerAngles(0, 1, 2);
    //   double pitch_data;
    //   pitch_data = euler[1];
    //   if(pitch_data>M_PI_2) pitch_data=M_PI-pitch_data;
    //   RCLCPP_INFO(this->get_logger(),"pitch:%f",pitch_data);

    //   // double pitch_data;
    //   // if(!this->joint_state_ptr_->position.empty()){
    //   //   pitch_data = this->joint_state_ptr_->position[0]; 
    //   // }

    
    //   double pitch_velo_data = this->imu_ptr_->angular_velocity.y; 
    //   pitch_pub.pitch = pitch_data+ 0.03*sin((imu_time-this->time_start)*10*2*M_PI); 
    //   pitch_pub.pitch_velo = pitch_velo_data+ 0.03/10*sin((imu_time-this->time_start)*10*2*M_PI); 
    //   RCLCPP_INFO(this->get_logger(),"pitch to be pub %f",pitch_pub.pitch);
    //   pitch_publisher->publish(pitch_pub);
    // }
    double pitch_data;
    double pitch_velo_data = this->imu_ptr_->angular_velocity.y; 
    if(!this->joint_state_ptr_->position.empty()){
      pitch_data = this->joint_state_ptr_->position[0]; 
      pitch_velo_data = this->imu_ptr_->angular_velocity.y; 
    }
    // pitch_pub.pitch = pitch_data; 
    // pitch_pub.pitch_velo = pitch_velo_data; 
    pitch_pub.pitch = pitch_data+ 0.05*sin((imu_time-this->time_start)*5*2*M_PI); 
    pitch_pub.pitch_velo = pitch_velo_data+ 0.05/5*sin((imu_time-this->time_start)*5*2*M_PI); 
    // RCLCPP_INFO(this->get_logger(),"pitch to be pub %f",pitch_pub.pitch);
    pitch_publisher->publish(pitch_pub);


    //3. extract x velocity from IMU to roll_velo
    roll_velo_pub.data = this->odom_ptr_->twist.twist.linear.x;
    roll_velo_publisher->publish(roll_velo_pub);

    //4. keep the contact state data and publish
    // this->contact_pub = this->contact;
    this->contact_pub.left_contact = this->contact_ptr_->left_contact;
    this->contact_pub.right_contact = this->contact_ptr_->right_contact;
    contact_publisher->publish(contact_pub);

    this->actuators_publihser->publish(this->actuator_cmds_pub);
    this->motion_publisher->publish(this->motion_cmds_pub);
    // RCLCPP_WARN(this->get_logger(),"get motion: %f",
    //     this->motion_cmds_pub.z_swing_des);
  }
  // callbacks for reading Simulator data
  void joint_state_callback(const sensor_msgs::msg::JointState::SharedPtr joint_state_msg)
  {
    this->joint_states.header.frame_id = joint_state_msg->header.frame_id;
    this->joint_states.header.stamp = joint_state_msg->header.stamp;

    // this->joint_states.position.resize(joint_state_msg->name.size());
    this->joint_state_ptr_->name.resize(joint_state_msg->name.size());
    this->joint_state_ptr_->position.resize(joint_state_msg->position.size());
    this->joint_state_ptr_->velocity.resize(joint_state_msg->velocity.size());
    for (size_t k = 0; k < this->joint_state_ptr_->name.size(); k++)
    {
      this->joint_state_ptr_->name[k] = joint_state_msg->name[k];
      this->joint_state_ptr_->position[k] = joint_state_msg->position[k];
      this->joint_state_ptr_->velocity[k] = joint_state_msg->velocity[k];
    }
    // for (size_t k =0; k<joint_state_msg->name.size(); k++){
    //     this->joint_states.position[k]=joint_state_msg->position[k];
    //     this->joint_states.velocity[k]=joint_state_msg->velocity[k];
    //     this->joint_states.name[k]=joint_state_msg->name[k];
    // }
    // RCLCPP_INFO(this->get_logger(),"%ld",this->joint_state_ptr_->position.size());
  }

  void imu_callback(const sensor_msgs::msg::Imu::SharedPtr imu_msg)
  {
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

  void contact_callback(const communication::msg::ContactState::SharedPtr contact_msg)
  {
    ////TODO: read in the contact states
    // this->contact = *contact_msg;
    this->contact_ptr_->header.frame_id = contact_msg->header.frame_id;
    this->contact_ptr_->header.stamp = contact_msg->header.stamp;

    this->contact_ptr_->left_contact.data = contact_msg->left_contact.data;
    this->contact_ptr_->right_contact.data = contact_msg->right_contact.data;

    // this->contact_pub.left_contact = contact_msg->left_contact.data;
    // this->contact_pub.right_contact = contact_msg->right_contact.data;
    
  }

  void odom_callback(const nav_msgs::msg::Odometry::SharedPtr odom_msg)
  {
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

  // void commands_callback()
  // {
  //   this->actuators_publihser->publish(this->actuator_cmds_pub);
  //   this->motion_publisher->publish(this->motion_cmds_pub);
  //   RCLCPP_WARN(this->get_logger(),"get motion: %f",
  //   this->motion_cmds_pub.z_swing_des);
  // }

  void actuators_callback(const communication::msg::JointState actuator_msgs){
    this->actuator_cmds_pub.actuators_name = {"right_hip","right_knee","left_hip","left_knee"};
    this->actuator_cmds_pub.torque_feedforward = {0,0,0,0};
    this->actuator_cmds_pub.kp = {this->kp_hip,this->kp_knee,this->kp_hip,this->kp_knee};
    this->actuator_cmds_pub.kd = {this->kd_hip,this->kd_knee,this->kd_hip,this->kd_knee};
    // this->actuator_cmds_pub.kp = {0,0,0,0};
    // this->actuator_cmds_pub.kd = {0,0,0,0};
    this->actuator_cmds_pub.pos_des.resize(4);
    this->actuator_cmds_pub.vel_des.resize(4);
    // this->actuator_cmds_pub.torque_feedforward.resize(4);
    for (size_t k =0; k<actuator_msgs.position.size(); k++){
      this->actuator_cmds_pub.torque_feedforward[k]=actuator_msgs.torque[k];
      this->actuator_cmds_pub.pos_des[k]=actuator_msgs.position[k];
      this->actuator_cmds_pub.vel_des[k]=actuator_msgs.velocity[k];
    }
    if(this->motion_cmds_pub.contact_foot == 0){ //left leg contact{
      this->actuator_cmds_pub.torque_feedforward[2]/=1.15;
      this->actuator_cmds_pub.torque_feedforward[3]/=1.15;
    }else{
      this->actuator_cmds_pub.torque_feedforward[0]/=0.75;
      this->actuator_cmds_pub.torque_feedforward[1]/=0.75;
    }
  }

  void motion_callback(const communication::msg::MotionCommands motion_msgs){
    this->motion_cmds_pub.com_height_des = motion_msgs.com_height_des;
    this->motion_cmds_pub.com_height_actual = motion_msgs.com_height_actual;
    this->motion_cmds_pub.x_des = motion_msgs.x_des;
    this->motion_cmds_pub.x_actual = motion_msgs.x_actual;
    this->motion_cmds_pub.z_swing_des = motion_msgs.z_swing_des;
    this->motion_cmds_pub.z_swing_actual = motion_msgs.z_swing_actual;
    this->motion_cmds_pub.pitch_des = motion_msgs.pitch_des;
    this->motion_cmds_pub.pitch_vel_des = motion_msgs.pitch_vel_des;
    this->motion_cmds_pub.pitch_vel_act = motion_msgs.pitch_vel_act;
    this->motion_cmds_pub.contact_foot = motion_msgs.contact_foot;
    this->motion_cmds_pub.vel_des.x= motion_msgs.vel_des.x;
    this->motion_cmds_pub.vel_des.y= motion_msgs.vel_des.y;
    this->motion_cmds_pub.vel_des.z= motion_msgs.vel_des.z;
  }

  void declare_params(){
    this->declare_parameter("kp_hip",10.);
    this->declare_parameter("kp_knee",10.);
    this->declare_parameter("kd_hip",0.1);
    this->declare_parameter("kd_knee",0.1);
    this->kp_hip = this->get_parameter("kp_hip").as_double();
    this->kp_knee = this->get_parameter("kp_knee").as_double();
    this->kd_hip = this->get_parameter("kd_hip").as_double();
    this->kd_knee = this->get_parameter("kd_knee").as_double();
    this->time_start = this->get_clock()->now().seconds();
  }

  rclcpp::TimerBase::SharedPtr timer_pub;
  rclcpp::TimerBase::SharedPtr actuators_cmds_timer;
  // declare publishers
  rclcpp::Publisher<communication::msg::PitchState>::SharedPtr pitch_publisher;
  rclcpp::Publisher<communication::msg::ContactState>::SharedPtr contact_publisher;
  rclcpp::Publisher<communication::msg::JointState>::SharedPtr joint_publisher;
  rclcpp::Publisher<std_msgs::msg::Float32>::SharedPtr roll_velo_publisher;

  rclcpp::Publisher<communication::msg::ActuatorCmds>::SharedPtr actuators_publihser;
  rclcpp::Publisher<communication::msg::MotionCommands>::SharedPtr motion_publisher;

  // declare subscribers
  rclcpp::Subscription<sensor_msgs::msg::JointState>::SharedPtr joint_state_subscription;
  rclcpp::Subscription<sensor_msgs::msg::Imu>::SharedPtr imu_subscription;
  rclcpp::Subscription<communication::msg::ContactState>::SharedPtr contact_subscription;
  rclcpp::Subscription<nav_msgs::msg::Odometry>::SharedPtr odom_subscription;

  rclcpp::Subscription<communication::msg::JointState>::SharedPtr actuators_subscription;
  rclcpp::Subscription<communication::msg::MotionCommands>::SharedPtr motion_subscription;  

  // size_t count_;
  std::string name_prefix;

  // message objects for storing subscribed data
  sensor_msgs::msg::JointState joint_states;
  sensor_msgs::msg::Imu IMU_data;
  communication::msg::ContactState contact;
  nav_msgs::msg::Odometry odom;
  // corresponding pointers for message objects
  std::shared_ptr<sensor_msgs::msg::JointState> joint_state_ptr_ = std::make_shared<sensor_msgs::msg::JointState>();
  std::shared_ptr<sensor_msgs::msg::Imu> imu_ptr_ = std::make_shared<sensor_msgs::msg::Imu>();
  std::shared_ptr<communication::msg::ContactState> contact_ptr_ = std::make_shared<communication::msg::ContactState>();
  std::shared_ptr<nav_msgs::msg::Odometry> odom_ptr_ = std::make_shared<nav_msgs::msg::Odometry>();

  std::shared_ptr<sensor_msgs::msg::JointState> actuators_cmds_sub_ptr_ = std::make_shared<sensor_msgs::msg::JointState>();
  std::shared_ptr<communication::msg::MotionCommands> motion_cmds_sub_ptr = std::make_shared<communication::msg::MotionCommands>();

  // message objects for storing data for publishing
  communication::msg::JointState joint_states_pub; 
  communication::msg::PitchState  pitch_pub; 
  std_msgs::msg::Float32  roll_velo_pub; 
  communication::msg::ContactState  contact_pub; 

  communication::msg::ActuatorCmds actuator_cmds_pub;
  communication::msg::MotionCommands motion_cmds_pub;

  float kp_hip,kp_knee,kd_hip,kd_knee;
  double time_start;
  ekf::EKF EKF_Filter;
};

int main(int argc, char *argv[])
{
  rclcpp::init(argc, argv);
  rclcpp::spin(std::make_shared<SimulatorToHardwarePublisher>());
  rclcpp::shutdown();
  return 0;
};