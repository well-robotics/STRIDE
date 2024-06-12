#include "rclcpp/rclcpp.hpp" 
#include "std_msgs/msg/float32.hpp"
#include "imu_lib.hpp"
#include "moving_average_filter.hpp"

#define BAR_LENGTH 1.30f

using namespace std::chrono_literals; 

/*************************This node is used to read IMU data through I2C bus****************************/
class RollVeloPublisher : public rclcpp::Node{
public: 
  RollVeloPublisher(): Node("roll_velo_publisher"),count_(0){
    imu1.sensor_address = 0x68;
    imu1.I2Cinit(1);

    
    auto qos = rclcpp::QoS(rclcpp::KeepLast(1), rmw_qos_profile_sensor_data);

    // initialize the publisher with a message type and topic name, queue size 
    publisher_ = this-> create_publisher<std_msgs::msg::Float32>("hardware/roll_velo",qos); 
    // setup the timer and the callback function 500ms measn one second run the callback function two times 
    timer_ = this -> create_wall_timer(5ms, std::bind(&RollVeloPublisher::timer_callback, this)); 
  }
private: 
  void timer_callback(){
    this->imu1.imuDataGet(); 
    this->MA_filter.renew((imu1.stGyroRawData.s16Z/GYRO_SCALE_FACTOR /180.0*M_PI )*BAR_LENGTH); 

    roll_velo.data = this->MA_filter.calc();//(imu1.stGyroRawData.s16X/GYRO_SCALE_FACTOR)*BAR_LENGTH;
    this->publisher_ -> publish(roll_velo); 
  }
rclcpp::TimerBase::SharedPtr timer_; 
rclcpp::Publisher<std_msgs::msg::Float32>::SharedPtr publisher_;
size_t count_; 
std_msgs::msg::Float32 roll_velo; 
// filter
MovingAverageFilter<double> MA_filter=MovingAverageFilter<double>(20,0.0); 
// the IMU device 
ICM20948 imu1; 
};
int main(int argc, char * argv[]){
  rclcpp::init(argc,argv); 
  rclcpp::spin(std::make_shared<RollVeloPublisher>()); 
  rclcpp::shutdown();
  return 0;
}