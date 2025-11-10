/*
 * @class: a ROS2 node for realizing H-LIP based walking
        It success the FiveLinkWalker class for planning the step 
        It sends desired output to lower level controller to follow the trajectory
 * @author: Yuhao Huang, Yicheng Zeng
*/

#include <S2S_controller.hpp>

struct Hlip_params{
    int orbit = 1; //default orbit=P1
    double vdes=1.0;
    double COM_height=0.3;
    double TSSP=0.5;
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

    float kp_hip = 20.0 ;
    float kp_knee = 20.0;
    float kd_hip = 0.1;
    float kd_knee = 0.1;
    double x_swing_test = 0.0;
    double z_swing_test = 0.0;
    double pitch_test = 0.0;
    double pelvis_init_height=0.450;
};

namespace hardware_ns
{
    class HlipController: public rclcpp::Node
    {
        public:
            HlipController(): Node("hlip_controller") ,name_prefix("hardware/")
            {
                RCLCPP_INFO(this->get_logger(), "controller node start, building publishers, subscribers and callbacks...");
                
                //declare parameters
                params_setup();

                //quality of service, buffer size of 1
                auto qos = rclcpp::QoS(rclcpp::KeepLast(1), rmw_qos_profile_sensor_data);

                // publishers to 
                // 1. joints states for hardware following
                // 2. motion commands for debugging
                actuators_cmds_publisher = this->create_publisher<communication::msg::JointState>(
                    name_prefix+"actuators_cmds", qos);
                motion_cmds_publisher = this->create_publisher<communication::msg::MotionCommands>(
                    name_prefix+"motion_cmds", qos);
                actuators_cmds_timer = this->create_wall_timer(
                    5ms, std::bind(&HlipController::actuators_cmds_callback,this));

                // subscribers to joint encoder, pitch data, COM_x velocity and contact sensSor
                joint_state_subscription = this->create_subscription<communication::msg::JointState>(
                    name_prefix+"joint_states", qos, std::bind(&HlipController::joint_state_callback, this, std::placeholders::_1));
                x_velo_subscription = this->create_subscription<std_msgs::msg::Float32>(
                    name_prefix+"roll_velo", qos, std::bind(&HlipController::x_velo_callback, this, std::placeholders::_1)); 
                imu_subscription = this->create_subscription<communication::msg::PitchState>(
                    name_prefix+"pitch",qos, std::bind( &HlipController::imu_callback, this, std::placeholders::_1 ));
                contact_subscription = this->create_subscription<communication::msg::ContactState>(
                    name_prefix+"ContactState", qos, std::bind( &HlipController::contact_callback, this, std::placeholders::_1 ));
                
                //timer_callback
                timer_gen = this->create_wall_timer(5ms, std::bind(&HlipController::timer_callback_hlip, this));

                pub_msg_setup();

                this->time_start = this->get_clock()->now().seconds();
                this->released = false;
                this->contact_ptr_->left_contact.data = false;
                this->contact_ptr_->right_contact.data = false;
                this->walker_class.output_controller_initialize(
                    this->hlip_params.orbit_init,this->hlip_params.vdes,this->hlip_params.COM_height,this->hlip_params.TSSP,
                    this->hlip_params.TDSP,this->hlip_params.step_size,this->hlip_params.swing_height,this->hlip_params.torso_angle);

            }
        private:

            void timer_callback_hlip(){
                this->time_now = this->get_clock()->now().seconds() -  this->time_start;
                this->actuator_cmds.position = {0.0, 0.0, 0.0, 0.0};
                this->actuator_cmds.velocity = {0.0, 0.0, 0.0, 0.0};

                // 1. read in the generalized coordinates q, where x and z are set to 
                VectorXd joint_pos = VectorXd::Zero(7);
                VectorXd joint_vel = VectorXd::Zero(7);
                VectorXd joint_state_fl = VectorXd::Zero(14);
                if(this->joint_state_ptr_->position.empty()){
                    // determine whether joints states are successfully read
                    // RCLCPP_WARN(this->get_logger(), "joint state a null pointer");
                    joint_pos << 0,0,0,0,0,0,0;
                    joint_vel << 0,0,0,0,0,0,0;
                    return;
                }else{
                    // 0,1,2: floating base x, z, pitch; 3,4,5,6: right_hip, right_knee, left_hip, left_knee
                    joint_pos(0) = 0;
                    joint_pos(1) = 0;
                    joint_pos(2) =  this->pitch;
                    for(int i=3; i<7; i++){ joint_pos(i)=this->joint_state_ptr_->position[i-3];}
                    // Reconfiguring the kinematics to dynamic library convention
                    joint_pos(3) = -joint_pos(3);
                    joint_pos(4) = M_PI - joint_pos(4);
                    joint_pos(5) = joint_pos(5);
                    joint_pos(6) = joint_pos(6)- M_PI;
                    //Velocity 
                    joint_vel(0) = this->x_velo_ptr_->data;
                    joint_vel(1) = 0;//this->x_velo_ptr_->data;
                    joint_vel(2) = this->pitch_vel ;
                    for(int i=3; i<7; i++){ joint_vel(i)= - this->joint_state_ptr_->velocity[i-3];}
                    joint_vel(5) = - joint_vel(5);
                    joint_vel(6) = - joint_vel(6);
                }
                joint_state_fl << joint_pos,joint_vel;

                // 2. determine if the robot started walking or not ){//
                if( this->time_now>=3.5  && !this->released   && (this->contact_ptr_->left_contact.data || this->contact_ptr_->right_contact.data)){ //
                    //release the robot once, when any of its feet first touch the ground
                    this->released=true;
                    RCLCPP_WARN(this->get_logger(),"Release the robot");
                }
                if( this->released && this->count_in_air>=200){
                    this->released=false;
                    this->count_in_air = 0;
                    this->actuator_cmds.position = {-0.6693, 0.9334, -0.6693, 0.9334 };
                    this->actuator_cmds.velocity = {0,0,0,0};
                    for(int i=0; i<4; i++){
                        this->actuator_cmds.torque[i] = 0.0;
                    }
                    RCLCPP_WARN(this->get_logger(),"back to initial configuration");
                    // sleep(0.5);
                    this->time_start=this->time_now;
                }
                // 2.1. If the robot was lifted from the ground for 1s, return to initial configuration.
                if( !this->contact_ptr_->left_contact.data && !this->contact_ptr_->right_contact.data){
                    this->count_in_air++;
                }else{
                    this->count_in_air=0;
                }

                // 2.2. Determine which leg should be the stance leg
                switch (this->contact_foot)
                {
                case left_leg:
                    if(
                        (this->time_now-this->walker_class.time_start >=  2.0 * this->walker_class.TSSP) || 
                        ((this->walker_class.T >= this->walker_class.TSSP )  && this->contact_ptr_->right_contact.data) ){//
                        this->previous_contact = left_leg;
                        this->contact_foot = right_leg;
                    }
                    break;
                case right_leg:
                    if(  
                        (this->time_now-this->walker_class.time_start >=  2.0 * this->walker_class.TSSP) || 
                        ((this->walker_class.T >= this->walker_class.TSSP)  && this->contact_ptr_->left_contact.data) ){//
                        this->previous_contact = right_leg;
                        this->contact_foot = left_leg;
                    }
                    break;
                default:
                    break;
                }

                // 3. hold the robot in zero configuration before feet touch the ground.
                if (this->time_now<2.0){    
                    // set the robot position to some initial configuration w.r.t. 0
                    // ensure IK does not fail
                    this->actuator_cmds.position = {-0.6693, 0.9334, -0.6693, 0.9334 };
                    this->actuator_cmds.velocity = {0,0,0,0};
                }
                else if(!this->released){
                    // Do pseudo renew of the hlip, just use ik to hold its configuration
                    Vector4d input_state_data = VectorXd::Zero(4); 
                    Vector4d d_input_state_data = VectorXd::Zero(4);
                    this->walker_class.SSP_output_des <<
                        this->hlip_params.COM_height,this->hlip_params.x_swing_test,
                        this->hlip_params.z_swing_test,this->hlip_params.pitch_test;
                    this->walker_class.SSP_output_vdes << 0,0,0,0;
                    // // The following code provides a sinesoidal wave test of endpoint
                    // // Useful to see if IK is going wrong.
                    // double phase = Pi/(0.5) * this->time_now;
                    // this->walker_class.SSP_output_des <<
                    //         this->hlip_params.COM_height,
                    //         0.05 * cos( phase ),
                    //         0.05 * (1-sin(phase)),
                    //         0.0;
                    // this->walker_class.SSP_output_vdes << 
                    //         0,
                    //         -0.05 * Pi/(0.5) * sin( phase ),
                    //         -0.05 * Pi/(0.5) * cos( phase ),
                    //         0;

                    MatrixXd COM_position_abs = MatrixXd::Zero(1,3);
                    COMPosition(COM_position_abs, joint_pos);
                    MatrixXd COMvelocity = MatrixXd::Zero(3,1);
                    COM_velocity(COMvelocity,joint_pos,joint_vel);
                    MatrixXd stance_leg_position = MatrixXd::Zero(2,1);
                    MatrixXd swing_leg_position_abs = VectorXd::Zero(2,1);
                    MatrixXd stance_leg_velocity = VectorXd::Zero(2,1);
                    MatrixXd swing_leg_velocity = VectorXd::Zero(2,1);
                    pRightToe(swing_leg_position_abs,joint_pos);
                    pLeftToe(stance_leg_position,joint_pos);
                    vRightToe(swing_leg_velocity, joint_state_fl);
                    vLeftToe(stance_leg_velocity,joint_state_fl);                    
                    input_state_data << (COM_position_abs(2) - stance_leg_position(1)),
                                        (swing_leg_position_abs(0)-stance_leg_position(0)),
                                        (swing_leg_position_abs(1)-stance_leg_position(1)),
                                        joint_pos(2);
                    d_input_state_data <<   COMvelocity(2) - stance_leg_velocity(1) ,
                                            swing_leg_velocity(0) - stance_leg_velocity(0) ,
                                            swing_leg_velocity(1) - stance_leg_velocity(1),
                                            joint_vel(2);

                    this->walker_class.own_q << joint_pos.segment<5>(2);
                    this->walker_class.output_controller_IK();
                    
                    VectorXd joint_pos_des = this->walker_class.output_controller_IK_solver.getSolution_Position();
                    VectorXd joint_vel_des = this->walker_class.output_controller_IK_solver.getSolution_Velocity();
                    
                    this->motion_cmds.com_height_des = this->walker_class.SSP_output_des(0); 
                    this->motion_cmds.com_height_actual = COM_position_abs(2) - stance_leg_position(1);// com_z - stance_z
                    this->motion_cmds.x_des = this->walker_class.SSP_output_des(1);
                    this->motion_cmds.x_actual = swing_leg_position_abs(0) - stance_leg_position(0);// swing_x - stance_x
                    this->motion_cmds.z_swing_des = this->walker_class.SSP_output_des(2);
                    this->motion_cmds.z_swing_actual = swing_leg_position_abs(1) - stance_leg_position(1);// swing_z - stance_z
                    this->motion_cmds.pitch_des = joint_pos(2);
                    this->motion_cmds.pitch_vel_des = joint_pos(2);
                    this->motion_cmds.pitch_vel_act = joint_vel(2);
                    this->motion_cmds.vel_des.x =this->walker_class.SSP_output_vdes(1);
                    this->motion_cmds.vel_des.z =this->walker_class.SSP_output_vdes(2);
                    this->motion_cmds.vel_des.y =this->walker_class.SSP_output_vdes(0);
                    this->motion_cmds.contact_foot = this->contact_foot==left_leg?0:1;

                    for(int i=0; i<4; i++){
                        this->actuator_cmds.position[i] = joint_pos_des(i+1);
                        this->actuator_cmds.velocity[i] = -joint_vel_des(i+1);  
                    }
                    // // match the configuration from dynamical lib to real robot
                    this->actuator_cmds.velocity[2] = -this->actuator_cmds.velocity[2];
                    this->actuator_cmds.velocity[3] = -this->actuator_cmds.velocity[3];
                    this->actuator_cmds.position[0] =    -    this->actuator_cmds.position[0];
                    this->actuator_cmds.position[1] = M_PI - this->actuator_cmds.position[1];
                    this->actuator_cmds.position[2] =    +  this->actuator_cmds.position[2];
                    this->actuator_cmds.position[3] = M_PI + this->actuator_cmds.position[3];
                }
                // 4. If the robot is detected to be released, it starts planning to walk
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
                    // COMvelocity(0)=joint_vel(0);
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

                    input_COM_state << (COM_position_abs(0) - stance_leg_position(0))// 
                                        ,COMvelocity(0) - stance_leg_velocity(0); //


                    //put in all the values calculated into HLIP
                    this->walker_class.output_controller_renew(joint_pos,joint_vel,
                    input_state_data,d_input_state_data,this->contact_foot,input_COM_state,this->time_now);

                    this->walker_class.output_controller_plan_desired();

                    // If due to any reason, the swing leg didn't touches the ground, make it streches down to find the ground.
                    if(this->time_now-this->walker_class.time_start > this->walker_class.TSSP){
                        this->walker_class.SSP_output_des(2)-= 0.2 * (this->time_now-this->walker_class.time_start - this->walker_class.TSSP) ;
                        this->walker_class.SSP_output_des(1)-= 2.0 * this->walker_class.SSP_output_des(1) *  (this->time_now-this->walker_class.time_start - this->walker_class.TSSP) ;
                    }

                    // send motion command for loggin and debugging
                    this->motion_cmds.com_height_des = this->walker_class.SSP_output_des(0);
                    this->motion_cmds.com_height_actual = COM_position_abs(2) - stance_leg_position(1);
                    this->motion_cmds.x_des = this->walker_class.SSP_output_des(1);
                    this->motion_cmds.x_actual = swing_leg_position_abs(0) - stance_leg_position(0);
                    this->motion_cmds.z_swing_des = this->walker_class.SSP_output_des(2);
                    this->motion_cmds.z_swing_actual = swing_leg_position_abs(1) - stance_leg_position(1);
                    this->motion_cmds.pitch_des = 0;
                    this->motion_cmds.pitch_act = joint_pos(2);
                    this->motion_cmds.com_x_act = (COM_position_abs(0) - stance_leg_position(0));
                    this->motion_cmds.com_x_des = (this->walker_class.COM_x_real(0));// - stance_leg_position(0));
                    this->motion_cmds.pitch_vel_des = joint_pos(2);
                    this->motion_cmds.pitch_vel_act = joint_vel(2);
                    this->motion_cmds.vel_des.x = joint_vel(0) - stance_leg_velocity(0);//this->walker_class.SSP_output_vdes(1);
                    this->motion_cmds.vel_des.z =this->walker_class.SSP_output_vdes(2);
                    this->motion_cmds.vel_des.y =this->walker_class.SSP_output_vdes(0);
                    this->motion_cmds.contact_foot = (this->contact_foot==left_leg)?0:1;

                    // Do IK to find desired joint position and velocity
                    this->walker_class.output_controller_IK();

                    if(this->walker_class.SSP_output_des.size() == 0){
                        RCLCPP_WARN(this->get_logger(), "Failed calculating IK output");
                    }
                    VectorXd joint_pos_des = this->walker_class.output_controller_IK_solver.getSolution_Position();
                    VectorXd joint_vel_des = this->walker_class.output_controller_IK_solver.getSolution_Velocity();
                    
                    // Do gravity compenstation 
                    VectorXd joint_pos_des_gravcomp = VectorXd::Zero(7);
                    joint_pos_des_gravcomp << 0., 0., joint_pos_des;
                    VectorXd tau_gravity = this->walker_class.output_controller_gravity_compensation(joint_pos_des_gravcomp);


                    for(int i=0; i<4; i++){
                        this->actuator_cmds.position[i] = joint_pos_des(i+1);
                        this->actuator_cmds.velocity[i] = -joint_vel_des(i+1); 
                        this->actuator_cmds.torque[i] = 0.80 *tau_gravity(i);
                    }
                    
                    // match the configuration from dynamical lib to real robot
                    this->actuator_cmds.velocity[2] = -this->actuator_cmds.velocity[2];
                    this->actuator_cmds.velocity[3] = -this->actuator_cmds.velocity[3];

                    this->actuator_cmds.torque[2] = -this->actuator_cmds.torque[2];
                    this->actuator_cmds.torque[3] = -this->actuator_cmds.torque[3];
                    if(this->contact_foot == left_leg){
                        this->actuator_cmds.torque[2] *= 1.10;
                        this->actuator_cmds.torque[3] *= 0.80;
                    }else{
                        this->actuator_cmds.torque[0] *= 0.75;
                        this->actuator_cmds.torque[1] *= 0.65;
                    }

                    for(int i=0; i<4; i++){
                        this->actuator_cmds.velocity[i]=0;
                    }
                    this->actuator_cmds.position[0] =    -    this->actuator_cmds.position[0];
                    this->actuator_cmds.position[1] = M_PI - this->actuator_cmds.position[1];
                    this->actuator_cmds.position[2] =    +  this->actuator_cmds.position[2];
                    this->actuator_cmds.position[3] = M_PI + this->actuator_cmds.position[3];
                    // RCLCPP_INFO(this->get_logger(),"step planned");
                    

                }
            }

            void actuators_cmds_callback(){ //publisher callback for actuators in form of ActuatorCmds msg
                this->actuators_cmds_publisher->publish(this->actuator_cmds);
                this->motion_cmds_publisher->publish(this->motion_cmds);
            }

            void joint_state_callback(const communication::msg::JointState::SharedPtr joint_state_msg){
                //read in list of joint position and velocity, in seq of [right_hip, right_knee, left_hip, left_knee]
                for (size_t k =0; k<4; k++){
                    // this->joint_state_ptr_->name[k]=joint_state_msg->name[k];
                    this->joint_state_ptr_->position[k]=joint_state_msg->position[k];
                    this->joint_state_ptr_->velocity[k]=joint_state_msg->velocity[k];
                }
            }

            void imu_callback(const communication::msg::PitchState pitch_msg){
                this->pitch = pitch_msg.pitch;
                this->pitch_vel = pitch_msg.pitch_velo;
                // this->pitch = 0;
                // this->pitch_vel = 0;
            }

            void contact_callback(const communication::msg::ContactState::SharedPtr contact_msg){
                //read in the contact state as 2 booleans
                this->contact_ptr_->header.frame_id = contact_msg->header.frame_id;
                this->contact_ptr_->header.stamp = contact_msg->header.stamp;

                this->contact_ptr_->left_contact.data = contact_msg->left_contact.data;
                this->contact_ptr_->right_contact.data = contact_msg->right_contact.data;
            }

            void x_velo_callback(const std_msgs::msg::Float32::SharedPtr x_velo_msg){
                this->x_velo_ptr_->data = x_velo_msg->data; 
            }

            void params_setup(){
                // some 
                this->declare_parameter("orbit_type",1); // 1 means P1 and 2 means P2
                this->declare_parameter("vdes",0.1); // vdes * Tssp = step_size for P1
                this->declare_parameter("TSSP",0.5);
                this->declare_parameter("TDSP",0.0);
                // the output variable of hlip
                this->declare_parameter("COM_height",0.30);
                this->declare_parameter("step_size",0.05); 
                this->declare_parameter("swing_height",0.015);
                this->declare_parameter("torso_angle", 0.);
                // initial z position of the floating base
                this->declare_parameter("pelvis_init_height",0.485);
                // motor pd parameters
                this->declare_parameter("kp_hip",10.);
                this->declare_parameter("kp_knee",10.);
                this->declare_parameter("kd_hip",0.1);
                this->declare_parameter("kd_knee",0.1);
                //initial IK members
                this->declare_parameter("x_swing_test",0.0);
                this->declare_parameter("z_swing_test",0.0);
                this->declare_parameter("pitch_test",0.0);
                // parameters for the hlip controller
                this->hlip_params.orbit = this->get_parameter("orbit_type").as_int();
                this->hlip_params.COM_height = this->get_parameter("COM_height").as_double();
                this->hlip_params.TSSP = this->get_parameter("TSSP").as_double();
                this->hlip_params.TDSP = this->get_parameter("TDSP").as_double();
                this->hlip_params.step_size = this->get_parameter("step_size").as_double();
                this->hlip_params.swing_height = this->get_parameter("swing_height").as_double();
                this->hlip_params.torso_angle = this->get_parameter("torso_angle").as_double();
                this->hlip_params.orbit_init = (this->hlip_params.orbit==1) ? P1 : P2;
                // I should reverse and let hlip determine step size itself, this is just for testing
                this->hlip_params.vdes = (this->hlip_params.orbit_init == P1)?
                                            this->hlip_params.step_size/this->hlip_params.TSSP :
                                            this->hlip_params.step_size/(2*this->hlip_params.TSSP);
                // init height
                this->hlip_params.pelvis_init_height = this->get_parameter("pelvis_init_height").as_double();
                // Motor pd parameters
                this->hlip_params.kp_hip = this->get_parameter("kp_hip").as_double();
                this->hlip_params.kp_knee = this->get_parameter("kp_knee").as_double();
                this->hlip_params.kd_hip = this->get_parameter("kd_hip").as_double();
                this->hlip_params.kd_knee = this->get_parameter("kd_knee").as_double();

                this->hlip_params.x_swing_test = this->get_parameter("x_swing_test").as_double();
                this->hlip_params.z_swing_test = this->get_parameter("z_swing_test").as_double();
                this->hlip_params.pitch_test = this->get_parameter("pitch_test").as_double();
            }

            void pub_msg_setup(){
                this->motion_cmds.com_height_actual=0.30;
                this->motion_cmds.com_height_des=this->hlip_params.COM_height;
                this->motion_cmds.x_actual=0.0;
                this->motion_cmds.x_des=0.0;
                this->motion_cmds.z_swing_actual=0.0;
                this->motion_cmds.z_swing_des=0.0;
                this->motion_cmds.pitch_des=0.0;
                this->motion_cmds.pitch_vel_act=0.0;
                this->motion_cmds.pitch_vel_des=0.0;
                this->motion_cmds.contact_foot=0;//left
                this->motion_cmds.vel_des.x=0.0;
                this->motion_cmds.vel_des.y=0.0;
                this->motion_cmds.vel_des.z=0.0;
            }

            void wrap(double &angle){
                //helper function to wrap up angles
                angle = fmod(angle+M_PI,2*M_PI)-M_PI;
            }

            Hlip_params hlip_params; // a struct wraps all the parameters
            output_controller walker_class;
            std::string name_prefix;
            stance_leg contact_foot=left_leg;
            stance_leg previous_contact = right_leg; //1 for left_foot, 2 for right
            double pitch=0;
            std::deque<double> pitch_list;
            double pitch_vel = 0;
            bool released=false;
            double time_now;
            double time_start;
            unsigned int count_in_air=0;

            // variable for IK smoothing
            VectorXd IK_filter_init = VectorXd::Zero(5,1); 
            MovingAverageFilter<VectorXd> IK_pos_filter = MovingAverageFilter<VectorXd>(15,IK_filter_init); 
            MovingAverageFilter<VectorXd> IK_velo_filter = MovingAverageFilter<VectorXd>(15,IK_filter_init);
            // MatrixXd joint_pos_des_prev = MatrixXd::Zero(5,2);
            // MatrixXd joint_vel_des_prev = MatrixXd::Zero(5,2);

            //the message sent
            communication::msg::JointState actuator_cmds;
            communication::msg::MotionCommands motion_cmds;
            
            //ptrs to be used in subscription callback
            std::shared_ptr<communication::msg::JointState> joint_state_ptr_  = std::make_shared<communication::msg::JointState>();
            std::shared_ptr<communication::msg::PitchState>        imu_ptr_          = std::make_shared<communication::msg::PitchState>();
            std::shared_ptr<communication::msg::ContactState> contact_ptr_  = std::make_shared<communication::msg::ContactState>();
            std::shared_ptr<std_msgs::msg::Float32> x_velo_ptr_ = std::make_shared<std_msgs::msg::Float32>(); 

            //shared pointers for timer callback.
            rclcpp::TimerBase::SharedPtr actuators_cmds_timer;
            rclcpp::TimerBase::SharedPtr timer_gen;
            //publishers
            rclcpp::Publisher<communication::msg::JointState>::SharedPtr actuators_cmds_publisher;
            rclcpp::Publisher<communication::msg::MotionCommands>::SharedPtr motion_cmds_publisher;
            //subscribers
            rclcpp::Subscription<communication::msg::JointState>::SharedPtr joint_state_subscription;
            rclcpp::Subscription<communication::msg::PitchState>::SharedPtr imu_subscription;
            rclcpp::Subscription<std_msgs::msg::Float32>::SharedPtr x_velo_subscription;
            rclcpp::Subscription<communication::msg::ContactState>::SharedPtr contact_subscription;
    };

} // namespace hardware_ns

int main(int argc, char * argv[])
{
    rclcpp::init(argc, argv);
    auto hlip_controller_node = std::make_shared<hardware_ns::HlipController>();
    rclcpp::spin(hlip_controller_node);
    rclcpp::shutdown();
    return 0;
}
