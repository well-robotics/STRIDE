#include "rclcpp/rclcpp.hpp" 
#include "sensor_msgs/msg/imu.hpp"
#include "std_msgs/msg/float32.hpp"
#include "imu_lib.hpp"
#include "madgwick_filter.hpp"
#include "communication/msg/pitch_state.hpp"
#include "moving_average_filter.hpp"

#include "low_pass_filter.hpp"

#undef I2C_ADD_ICM20948
#define I2C_ADD_ICM20948 0x69

using namespace std::chrono_literals; 

/*************************This node is used to read IMU data through I2C bus****************************/
class IMUPublisher : public rclcpp::Node{
public: 
  IMUPublisher(): Node("imu_publisher"),count_(0){
    imu1.sensor_address = 0x69; 
    imu1.I2Cinit(1);

    auto qos = rclcpp::QoS(rclcpp::KeepLast(1), rmw_qos_profile_sensor_data);

    // initialize the publisher with a message type and topic name, queue size 
    // publisher_ = this -> create_publisher<sensor_msgs::msg::Imu>("imu/data_raw",10); 
    //pitch_publisher = this -> create_publisher<std_msgs::msg::Float32>("hardware/pitch",qos); 
    pitch_publisher = this -> create_publisher<communication::msg::PitchState>("hardware/pitch",qos); 
    // setup the timer and the callback function 500ms measn one second run the callback function two times 
    timer_ = this -> create_wall_timer(5ms, std::bind(&IMUPublisher::timer_callback, this)); 
    exp_filter.set_coeff(0.9);

    this->time_start = this->get_clock()->now().seconds();
  }
private: 
  void timer_callback(){
    this->imu1.imuDataGet(); 

    /*********************directly use madgwick filter */ 
    // TODO: publish both the pitch angle and pitch velocity
    double ax = imu1.stAccelRawData.s16X/ACCEL_SCALE_FACTOR;
    double ay = imu1.stAccelRawData.s16Y/ACCEL_SCALE_FACTOR;
    double az = imu1.stAccelRawData.s16Z/ACCEL_SCALE_FACTOR;
    //units of gyroscope is in degrees/s, switch to rad/s first
    double gx = imu1.stGyroRawData.s16X/GYRO_SCALE_FACTOR *M_PI/180.; 
    double gy = imu1.stGyroRawData.s16Y/GYRO_SCALE_FACTOR *M_PI/180.;
    double gz = imu1.stGyroRawData.s16Z/GYRO_SCALE_FACTOR *M_PI/180.;

    //************Using EKF filter to filter orientation************//
    // double imu_time = static_cast<double>(rclcpp::Clock().now().nanoseconds()) / 1e9 - this->time_start;
    // Vector3d accel = Vector3d::Zero();
    // Vector3d angular_vel = Vector3d::Zero();
    // accel << ax, ay, az;
    // angular_vel << gx, gy, gz;
    // this->EKF_filter.update_raw_data(accel,angular_vel,imu_time);
    // Quaternionf q = this->EKF_filter.predict(angular_vel);
    // q = EKF_filter.correct(accel);
    // auto euler = q.toRotationMatrix().eulerAngles(0, 1, 2);

    // // Quaternionf q(imu1.stQuaternion.q0, imu1.stQuaternion.q1, imu1.stQuaternion.q2, imu1.stQuaternion.q3);
    // // auto euler = q.toRotationMatrix().eulerAngles(0, 1, 2);
    
    // this->pitch.pitch = euler(1);
    // if(this->pitch.pitch > M_PI_2) this->pitch.pitch = M_PI - this->pitch.pitch;
    // if(this->pitch.pitch < -M_PI_2) this->pitch.pitch = -this->pitch.pitch - M_PI;
    // this->pitch.pitch_velo = gy;
    //************end************//

    // //************Using madgewick filter to filter orientation************//
    this->filter.set_data(ax,ay,az,gx,gy,gz); 
    this->filter.update(); 
    this->filter.quaternion_to_Euler();
    
    this->pitch.pitch = this->filter.get_pitch_angle();
    this->pitch.pitch_velo = gy;
    // //************end************//

    // smooth the output with moving average filter
    // renew
    // this->average_filter[0].renew(this->filter.get_pitch_angle()); 
    // this->average_filter[1].renew(gy); 
    // // calc
    // this->pitch.pitch = this->average_filter[0].calc(); 
    // this->pitch.pitch_velo = this->average_filter[1].calc(); 


    this->pitch.pitch = exp_filter.calc(this->pitch.pitch);

    this->pitch_publisher -> publish(this->pitch);
  }

  rclcpp::TimerBase::SharedPtr timer_; 
  rclcpp::Publisher<communication::msg::PitchState>::SharedPtr pitch_publisher;
  size_t count_; 

  // pitch data
  // sensor_msgs::msg::Imu imu_raw; 
  // std_msgs::msg::Float32 pitch; 
  communication::msg::PitchState pitch;

  // the IMU reading object
  ICM20948 imu1;
  // the filter object 
  double dt = 0.001; 
  double gyro_e = 0; 
  double b = 0;
  Madgwick_filter filter; 
  MovingAverageFilter<double> average_filter[2] = {MovingAverageFilter<double>(10,0.0), MovingAverageFilter<double>(10,0.0)}; // average_filter[0] for pitch angle and average_filter[1] for pitch velocity
  LowPassFilter exp_filter;

  double time_start;

};
int main(int argc, char * argv[]){
  rclcpp::init(argc,argv); 
  rclcpp::spin(std::make_shared<IMUPublisher>()); 
  rclcpp::shutdown();
  return 0;
}