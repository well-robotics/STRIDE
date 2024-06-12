// custom code for ROS2

#include "MuJoCoMessageHandler.h"

namespace mujoco_ros2
{

  MuJoCoMessageHandler::MuJoCoMessageHandler(mj::Simulate *sim)
      : Node("MuJoCoMessageHandler"), sim_(sim), name_prefix("simulation/")
  {
    // set up name
    // model_param_name = name_prefix + "model_file";
    model_param_name = "model_file";

    this->declare_parameter(model_param_name, "");

    // set up service
    reset_service_ = this->create_service<communication::srv::SimulationReset>(
        name_prefix + "sim_reset",
        std::bind(&MuJoCoMessageHandler::reset_callback, this,
                  std::placeholders::_1, std::placeholders::_2));

    // communication quality set up
    auto qos = rclcpp::QoS(rclcpp::KeepLast(1), rmw_qos_profile_sensor_data);

    imu_publisher_ = this->create_publisher<sensor_msgs::msg::Imu>(
        name_prefix + "imu_data", qos);
    joint_state_publisher_ = this->create_publisher<sensor_msgs::msg::JointState>(
        name_prefix + "joint_states", qos);
    actuator_state_publisher_ = this->create_publisher<sensor_msgs::msg::JointState>(
        name_prefix + "actuator_states", qos);
    contact_sensor_publisher_ = this->create_publisher<communication::msg::ContactForce>(
        name_prefix + "contact_states", qos);
    odom_publisher_ = this->create_publisher<nav_msgs::msg::Odometry>(
        name_prefix + "odom", qos);
    sim_time_publisher_ = this->create_publisher<std_msgs::msg::Float64>(
        name_prefix + "sim_time", qos);
    //Added Feb 06
    contact_state_publisher_= this->create_publisher<communication::msg::ContactState>(
        name_prefix + "ContactState", qos);

    //Added April 02
    // ctrl_out_publisher_ = this->create_publisher<communication::msg::JointState>(
    //     name_prefix+ "sim_ctrl",qos);

    // Frequency setup by create a timer with callback function
    // IMU 400hz
    timers_.emplace_back(this->create_wall_timer(
        2.5ms, std::bind(&MuJoCoMessageHandler::imu_callback, this)));
    // joint states 1000hz; for debugging purposes, not provided on physical robot
    timers_.emplace_back(this->create_wall_timer(
        1ms, std::bind(&MuJoCoMessageHandler::joint_callback, this)));
    // encoder 1000hz
    timers_.emplace_back(this->create_wall_timer(
        1ms, std::bind(&MuJoCoMessageHandler::actuator_feedback_callback, this)));


    timers_.emplace_back(this->create_wall_timer(
        1ms, std::bind(&MuJoCoMessageHandler::contact_callback, this)));

    // time sim 1000hz
    timers_.emplace_back(this->create_wall_timer(
        1ms, std::bind(&MuJoCoMessageHandler::sim_time_callback, this)));

    // contact 400hz
    timers_.emplace_back(this->create_wall_timer(
        2.5ms, std::bind(&MuJoCoMessageHandler::contact_state_callback, this)));

    // odom 50hz or 400HZ, really depending on what you do
    timers_.emplace_back(this->create_wall_timer(
        2.5ms, std::bind(&MuJoCoMessageHandler::odom_callback, this)));
    // drop old 10hz
    timers_.emplace_back(this->create_wall_timer(
        100ms, std::bind(&MuJoCoMessageHandler::drop_old_message, this)));

    // subscribe actuator command
    actuator_cmd_subscription_ =
        this->create_subscription<communication::msg::ActuatorCmds>(
            // "bucky/actuators_cmds", qos,
            "simulation/actuators_cmds", qos,
            std::bind(&MuJoCoMessageHandler::actuator_cmd_callback, this,
                      std::placeholders::_1));

    param_subscriber_ = std::make_shared<rclcpp::ParameterEventHandler>(this);
    cb_handle_ = param_subscriber_->add_parameter_callback(
        model_param_name, std::bind(&MuJoCoMessageHandler::parameter_callback,
                                    this, std::placeholders::_1));

    actuator_cmds_ptr_ = std::make_shared<ActuatorCmds>();

    // RCLCPP_INFO(this->get_logger(), "Start MuJoCoMessageHandler ...");

    std::string model_file = this->get_parameter(model_param_name)
                                 .get_parameter_value()
                                 .get<std::string>();
    mju::strcpy_arr(sim_->filename, model_file.c_str());

    //
    sim_->uiloadrequest.fetch_add(1);

    //
    // this->init_logging();
  }

  MuJoCoMessageHandler::~MuJoCoMessageHandler()
  {
    // RCLCPP_INFO(this->get_logger(), "close node ...");
  }

  void MuJoCoMessageHandler::reset_callback(
      const std::shared_ptr<communication::srv::SimulationReset::Request> request,
      std::shared_ptr<communication::srv::SimulationReset::Response> response)
  {
    while (sim_->d == nullptr && rclcpp::ok())
    {
      std::this_thread::sleep_for(std::chrono::seconds(1));
    }
    if (sim_->d != nullptr)
    {
      if (request->header.frame_id != std::string(&sim_->m->names[0]))
      {
        // RCLCPP_ERROR(this->get_logger(), "reset request is not for %s",
        //              &sim_->m->names[0]);
        response->is_success = false;
      }
      else
      {
        sim_->mtx.lock();
        mj_resetData(sim_->m, sim_->d);
        sim_->d->qpos[0] = request->base_pose.position.x;
        sim_->d->qpos[1] = request->base_pose.position.y;
        sim_->d->qpos[2] = request->base_pose.position.z;
        sim_->d->qpos[3] = request->base_pose.orientation.w;
        sim_->d->qpos[4] = request->base_pose.orientation.x;
        sim_->d->qpos[5] = request->base_pose.orientation.y;
        sim_->d->qpos[6] = request->base_pose.orientation.z;

        for (int i = 0; i < request->joint_state.position.size(); i++)
        {
          int joint_id = mj_name2id(sim_->m, mjOBJ_JOINT,
                                    request->joint_state.name[i].c_str());
          if (joint_id > -1)
          {
            sim_->d->qpos[sim_->m->jnt_qposadr[joint_id]] =
                request->joint_state.position[i];
          }
          else
          {
            // RCLCPP_WARN(this->get_logger(), "Request joint %s does not exist",
            //             request->joint_state.name[i].c_str());
          }
        }
        for (size_t k = 0; k < actuator_cmds_ptr_->actuators_name.size(); k++)
        {
          actuator_cmds_ptr_->kp[k] = 0;
          actuator_cmds_ptr_->pos_des[k] = 0;
          actuator_cmds_ptr_->kd[k] = 0;
          actuator_cmds_ptr_->vel_des[k] = 0;
          actuator_cmds_ptr_->torque_feedforward[k] = 0.0;
        }
        sim_->mtx.unlock();
        response->is_success = true;
        // RCLCPP_INFO(this->get_logger(), "reset robot state...");
        // RCLCPP_INFO(this->get_logger(), "robot total mass: %f",
                    // sim_->m->body_subtreemass[0]);
      }
    }
    else
    {
      response->is_success = false;
    }
  }

  void MuJoCoMessageHandler::imu_callback()
  {
    if (sim_->d != nullptr)
    {
      auto message = sensor_msgs::msg::Imu();
      message.header.frame_id = &sim_->m->names[0];
      message.header.stamp = rclcpp::Clock().now();
      const std::lock_guard<std::mutex> lock(sim_->mtx);

      std::random_device rd;
      std::mt19937 gen(rd());
      std::normal_distribution<double> dist(0.0,1.0);

      for (int i = 0; i < sim_->m->nsensor; i++)
      {
        if (sim_->m->sensor_type[i] == mjtSensor::mjSENS_ACCELEROMETER)
        {
          message.linear_acceleration.x =
              sim_->d->sensordata[sim_->m->sensor_adr[i]];
          message.linear_acceleration.y =
              sim_->d->sensordata[sim_->m->sensor_adr[i] + 1];
          message.linear_acceleration.z =
              sim_->d->sensordata[sim_->m->sensor_adr[i] + 2];
        }
        else if (sim_->m->sensor_type[i] == mjtSensor::mjSENS_FRAMEQUAT)
        {
          message.orientation.w = sim_->d->sensordata[sim_->m->sensor_adr[i]];
          message.orientation.x = sim_->d->sensordata[sim_->m->sensor_adr[i] + 1];
          message.orientation.y = sim_->d->sensordata[sim_->m->sensor_adr[i] + 2];
          message.orientation.z = sim_->d->sensordata[sim_->m->sensor_adr[i] + 3];
        }
        else if (sim_->m->sensor_type[i] == mjtSensor::mjSENS_GYRO)
        {
          message.angular_velocity.x =
              sim_->d->sensordata[sim_->m->sensor_adr[i]];
          message.angular_velocity.y =
              sim_->d->sensordata[sim_->m->sensor_adr[i] + 1];
          message.angular_velocity.z =
              sim_->d->sensordata[sim_->m->sensor_adr[i] + 2];
        }
      }
      imu_publisher_->publish(message);
    }
  }

  void MuJoCoMessageHandler::odom_callback()
  {
    const std::lock_guard<std::mutex> lock(sim_->mtx);
    if (sim_->d != nullptr)
    {
      auto message = nav_msgs::msg::Odometry();
      message.header.frame_id = &sim_->m->names[0];
      message.header.stamp = rclcpp::Clock().now();
      message.pose.pose.position.x = sim_->d->qpos[0];
      message.pose.pose.position.y = sim_->d->qpos[1];
      message.pose.pose.position.z = sim_->d->qpos[2];
      message.pose.pose.orientation.w = sim_->d->qpos[3];
      message.pose.pose.orientation.x = sim_->d->qpos[4];
      message.pose.pose.orientation.y = sim_->d->qpos[5];
      message.pose.pose.orientation.z = sim_->d->qpos[6];
      message.twist.twist.linear.x = sim_->d->qvel[0];
      message.twist.twist.linear.y = sim_->d->qvel[1];
      message.twist.twist.linear.z = sim_->d->qvel[2];
      message.twist.twist.angular.x = sim_->d->qvel[3];
      message.twist.twist.angular.y = sim_->d->qvel[4];
      message.twist.twist.angular.z = sim_->d->qvel[5];
      odom_publisher_->publish(message);
    }
  }

  void MuJoCoMessageHandler::joint_callback()
  { //  publish all joint states.
    const std::lock_guard<std::mutex> lock(sim_->mtx);
    if (sim_->d != nullptr)
    {
      sensor_msgs::msg::JointState jointState;
      jointState.header.frame_id = &sim_->m->names[0];
      jointState.header.stamp = rclcpp::Clock().now();
      for (int i = 0; i < sim_->m->njnt; i++)
      {
        if (sim_->m->jnt_type[i] == mjtJoint::mjJNT_HINGE)
        {
          std::string jnt_name(mj_id2name(sim_->m, mjtObj::mjOBJ_JOINT, i));
          jointState.name.emplace_back(jnt_name);
          jointState.position.push_back(sim_->d->qpos[sim_->m->jnt_qposadr[i]]);
          jointState.velocity.push_back(sim_->d->qvel[sim_->m->jnt_dofadr[i]]);
          jointState.effort.push_back(
              sim_->d->qfrc_actuator[sim_->m->jnt_dofadr[i]]);

          // std::cout << "joint name: " << jnt_name << "  "
          //           << sim_->d->qpos[sim_->m->jnt_qposadr[i]]
          //           << std::endl;
        }
      }
      joint_state_publisher_->publish(jointState);
    }
  }

  void MuJoCoMessageHandler::sim_time_callback()
  {
    const std::lock_guard<std::mutex> lock(sim_->mtx);
    if (sim_->d != nullptr)
    {
      std_msgs::msg::Float64 msg;
      msg.data = sim_->d->time;
      simTime = sim_->d->time;
      sim_time_publisher_->publish(msg);
    }
    // logger.spin_logging();
  }

  void MuJoCoMessageHandler::actuator_feedback_callback()
  { // publish joint with encoders, not actual sensor value unless all sensors are put in bucky.xml
    const std::lock_guard<std::mutex> lock(sim_->mtx);
    if (sim_->d != nullptr)
    {
      sensor_msgs::msg::JointState actuatorState;
      actuatorState.header.frame_id = &sim_->m->names[0]; // robot name
      actuatorState.header.stamp = rclcpp::Clock().now();
      /* for (int i = 0; i < sim_->m->nsensor; i++)
       { // use actuatorpos to locate its joint, extract the sensor value,
         // sensor joint name -> id-> velocity and actuator force
         // note actuator (joint sensing) force cannot be used unless the robot has joint ft sensor
         if (sim_->m->sensor_type[i] == mjSENS_ACTUATORPOS)
         {
           std::string actuator_name(mj_id2name(sim_->m, mjtObj::mjOBJ_ACTUATOR, sim_->m->sensor_objid[i]));
           actuator_name.erase(0, 5); // by the convection, remove first 5 get the joint name;
           int joint_id(mj_name2id(sim_->m, mjtObj::mjOBJ_JOINT, actuator_name.c_str()));
           encoderState.name.emplace_back(actuator_name);
           encoderState.position.emplace_back(sim_->d->qpos[sim_->m->jnt_qposadr[joint_id]]);
           encoderState.velocity.emplace_back(sim_->d->qvel[sim_->m->jnt_dofadr[joint_id]]);
           encoderState.effort.emplace_back(
               sim_->d->qfrc_actuator[sim_->m->jnt_dofadr[joint_id]]);
         }
       }
       encoder_state_publisher_->publish(encoderState);
     }*/
      // ALTERNATIVE: define actuatorpos, actuatorvel, actuatorfrc
      // use string operation, which requires convections in
      for (int i = 0; i < sim_->m->nsensor; i++)
      {
        if (sim_->m->sensor_type[i] == mjSENS_ACTUATORPOS)
        {

          // s.t. right_hip_roll_joint_pos
          std::string actuator_pos_name = mj_id2name(sim_->m, mjtObj::mjOBJ_SENSOR, sim_->m->sensor_objid[i]);
          int pos_sensor_id = mj_name2id(sim_->m, mjOBJ_SENSOR,
                                         (actuator_pos_name).c_str());
          // s.t. right_hip_roll_joint
          std::string actuator_name = actuator_pos_name.substr(0, actuator_pos_name.length() - 4);
          int vel_sensor_id = mj_name2id(sim_->m, mjOBJ_SENSOR,
                                         (actuator_name + "_vel").c_str());
          int frc_sensor_id = mj_name2id(sim_->m, mjOBJ_SENSOR,
                                         (actuator_name + "_frc").c_str());
          actuatorState.name.emplace_back(actuator_name);
          actuatorState.position.emplace_back(sim_->d->sensordata[sim_->m->sensor_adr[pos_sensor_id]]);
          actuatorState.velocity.emplace_back(sim_->d->sensordata[sim_->m->sensor_adr[vel_sensor_id]]);
          actuatorState.effort.emplace_back(sim_->d->sensordata[sim_->m->sensor_adr[frc_sensor_id]]);
        }
      }

      actuator_state_publisher_->publish(actuatorState);
    }
  }


  void MuJoCoMessageHandler::contact_state_callback()
  {
    //publish contact state of both legs
    const std::lock_guard<std::mutex> lock(sim_->mtx);
    if (sim_->d != nullptr){
      communication::msg::ContactState contactstate;
      double threshold = 2.0;
      contactstate.header.frame_id = &sim_->m->names[0];//robot's name
      contactstate.header.stamp = rclcpp::Clock().now();
      double temp_force[1];
      std::string sensor_name = "left_contact";
      mju_copy(temp_force, &sim_->d->sensordata[sim_->m->sensor_adr[mj_name2id(sim_->m, mjOBJ_SENSOR, sensor_name.c_str())]], 1);

      contactstate.left_contact.data = (temp_force[0]>threshold) ? true : false;

      sensor_name = "right_contact";
      mju_copy(temp_force, &sim_->d->sensordata[sim_->m->sensor_adr[mj_name2id(sim_->m, mjOBJ_SENSOR, sensor_name.c_str())]], 1);
      contactstate.right_contact.data = (temp_force[0]>threshold) ? true : false;
      // RCLCPP_INFO(rclcpp::get_logger("MuJoCo"),"right contact forces: %f",temp_force[0]);
      contact_state_publisher_->publish(contactstate);
    }
  }

  void MuJoCoMessageHandler::contact_callback()
  { // publish joint with encoders, not actual sensor value unless all sensors are put in bucky.xml
    const std::lock_guard<std::mutex> lock(sim_->mtx);
    if (sim_->d != nullptr)
    {

      communication::msg::ContactForce contactForce;
      contactForce.header.frame_id = &sim_->m->names[0]; // robot name
      contactForce.header.stamp = rclcpp::Clock().now();
      // use string operation, which requires convections in xml
      std::string force_sensor_name = "right_inner_toe";
      double temp_forces[3];
      int sensor_id = mj_name2id(sim_->m, mjOBJ_SENSOR, force_sensor_name.c_str());

      mju_copy(temp_forces, &sim_->d->sensordata[sim_->m->sensor_adr[mj_name2id(sim_->m, mjOBJ_SENSOR, force_sensor_name.c_str())]], 3);

      force_sensor_name = "right_inner_heel";
      mju_copy(temp_forces, &sim_->d->sensordata[sim_->m->sensor_adr[mj_name2id(sim_->m, mjOBJ_SENSOR, force_sensor_name.c_str())]], 3);
      // force sensor is the sandwich between the two bodies on the robot.
      // Forces are on the opposite direction.
      contactForce.right_inner_heel.x = -temp_forces[0];
      contactForce.right_inner_heel.y = -temp_forces[1];
      contactForce.right_inner_heel.z = -temp_forces[2];

      force_sensor_name = "right_outer_toe";
      mju_copy(temp_forces, &sim_->d->sensordata[sim_->m->sensor_adr[mj_name2id(sim_->m, mjOBJ_SENSOR, force_sensor_name.c_str())]], 3);
      contactForce.right_outer_toe.x = -temp_forces[0];
      contactForce.right_outer_toe.y = -temp_forces[1];
      contactForce.right_outer_toe.z = -temp_forces[2];

      force_sensor_name = "right_inner_toe";
      mju_copy(temp_forces, &sim_->d->sensordata[sim_->m->sensor_adr[mj_name2id(sim_->m, mjOBJ_SENSOR, force_sensor_name.c_str())]], 3);
      contactForce.right_inner_toe.x = -temp_forces[0];
      contactForce.right_inner_toe.y = -temp_forces[1];
      contactForce.right_inner_toe.z = -temp_forces[2];

      force_sensor_name = "right_outer_heel";
      mju_copy(temp_forces, &sim_->d->sensordata[sim_->m->sensor_adr[mj_name2id(sim_->m, mjOBJ_SENSOR, force_sensor_name.c_str())]], 3);
      contactForce.right_outer_heel.x = -temp_forces[0];
      contactForce.right_outer_heel.y = -temp_forces[1];
      contactForce.right_outer_heel.z = -temp_forces[2];

      force_sensor_name = "left_outer_toe";
      mju_copy(temp_forces, &sim_->d->sensordata[sim_->m->sensor_adr[mj_name2id(sim_->m, mjOBJ_SENSOR, force_sensor_name.c_str())]], 3);
      contactForce.left_outer_toe.x = -temp_forces[0];
      contactForce.left_outer_toe.y = -temp_forces[1];
      contactForce.left_outer_toe.z = -temp_forces[2];

      force_sensor_name = "left_outer_heel";
      mju_copy(temp_forces, &sim_->d->sensordata[sim_->m->sensor_adr[mj_name2id(sim_->m, mjOBJ_SENSOR, force_sensor_name.c_str())]], 3);
      contactForce.left_outer_heel.x = -temp_forces[0];
      contactForce.left_outer_heel.y = -temp_forces[1];
      contactForce.left_outer_heel.z = -temp_forces[2];

      force_sensor_name = "left_inner_toe";
      mju_copy(temp_forces, &sim_->d->sensordata[sim_->m->sensor_adr[mj_name2id(sim_->m, mjOBJ_SENSOR, force_sensor_name.c_str())]], 3);
      contactForce.left_inner_toe.x = -temp_forces[0];
      contactForce.left_inner_toe.y = -temp_forces[1];
      contactForce.left_inner_toe.z = -temp_forces[2];

      force_sensor_name = "left_inner_heel";
      mju_copy(temp_forces, &sim_->d->sensordata[sim_->m->sensor_adr[mj_name2id(sim_->m, mjOBJ_SENSOR, force_sensor_name.c_str())]], 3);
      contactForce.left_inner_heel.x = -temp_forces[0];
      contactForce.left_inner_heel.y = -temp_forces[1];
      contactForce.left_inner_heel.z = -temp_forces[2];
      contact_sensor_publisher_->publish(contactForce);
    }
  }

  void MuJoCoMessageHandler::actuator_cmd_callback(
      const communication::msg::ActuatorCmds::SharedPtr msg) const
  {

    if (sim_->d != nullptr)
    {
      actuator_cmds_ptr_->time = (this->now().nanoseconds() - time_start) * 1e-9;

      // std::cout << "ros time:" << actuator_cmds_ptr_->time << std::endl;
      // parse control message
      // RCLCPP_INFO(rclcpp::get_logger("MuJoCo"), "Receive actuator commands, %ld %s", msg->actuators_name.size(), msg->actuators_name[0].c_str() );
      actuator_cmds_ptr_->actuators_name.resize(msg->actuators_name.size());
      actuator_cmds_ptr_->kp.resize(msg->kp.size());
      actuator_cmds_ptr_->pos_des.resize(msg->pos_des.size());
      actuator_cmds_ptr_->kd.resize(msg->kd.size());
      actuator_cmds_ptr_->vel_des.resize(msg->vel_des.size());
      actuator_cmds_ptr_->torque_feedforward.resize(msg->torque_feedforward.size());
      // RCLCPP_INFO(rclcpp::get_logger("MuJoCo"), "%ld %ld %ld %ld %ld %ld",
      //   actuator_cmds_ptr_->actuators_name.size(),actuator_cmds_ptr_->kp.size(),actuator_cmds_ptr_->pos_des.size(),
      //   actuator_cmds_ptr_->kd.size(),actuator_cmds_ptr_->vel_des.size(),actuator_cmds_ptr_->torque_feedforward.size()  );
      for (size_t k = 0; k < msg->actuators_name.size(); k++)
      {
        actuator_cmds_ptr_->actuators_name[k] = msg->actuators_name[k];
        actuator_cmds_ptr_->kp[k] = msg->kp[k];
        actuator_cmds_ptr_->pos_des[k] = msg->pos_des[k];
        actuator_cmds_ptr_->kd[k] = msg->kd[k];
        actuator_cmds_ptr_->vel_des[k] = msg->vel_des[k];
        actuator_cmds_ptr_->torque_feedforward[k] = msg->torque_feedforward[k];
      }


      // RCLCPP_INFO(this->get_logger(), "subscribe actuator cmds");
    }
  }

  void MuJoCoMessageHandler::parameter_callback(const rclcpp::Parameter &)
  {
    std::string model_file = this->get_parameter(model_param_name)
                                 .get_parameter_value()
                                 .get<std::string>();
    // RCLCPP_INFO(this->get_logger(), "load model from: %s", model_file.c_str());
    mju::strcpy_arr(sim_->filename, model_file.c_str());
    sim_->uiloadrequest.fetch_add(1);
  }

  void MuJoCoMessageHandler::drop_old_message()
  {
    // check if it has been a while since last received command
    if (abs(actuator_cmds_ptr_->time - this->now().seconds()) > 0.2)
    {
      for (size_t k = 0; k < actuator_cmds_ptr_->actuators_name.size(); k++)
      {
        actuator_cmds_ptr_->kp[k] = 0.0;
        actuator_cmds_ptr_->pos_des[k] = 0.0;
        actuator_cmds_ptr_->kd[k] = 1.0;
        actuator_cmds_ptr_->vel_des[k] = 0.0;
        actuator_cmds_ptr_->torque_feedforward[k] = 0.0;
      }
    }
  }

  std::shared_ptr<MuJoCoMessageHandler::ActuatorCmds>
  MuJoCoMessageHandler::get_actuator_cmds_ptr()
  {
    return actuator_cmds_ptr_;
  }

} // namespace mujoco_ros2
