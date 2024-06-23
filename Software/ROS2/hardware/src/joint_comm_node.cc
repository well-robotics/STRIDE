#include <memory>
#include <functional>
#include <libserial/SerialPort.h>
#include <libserial/SerialStream.h>
#include <iostream>
#include <chrono>
#include <thread>
#include <unistd.h>
#include <stdint.h>

#include "rclcpp/rclcpp.hpp"
//#include "sensor_msgs/msg/joint_state.hpp"
#include "communication/msg/joint_state.hpp"
#include "low_pass_filter.hpp"

using namespace LibSerial;
using namespace std;
using std::placeholders::_1;
using namespace std::chrono_literals;

class SerialJointNode : public rclcpp::Node
{
public:
  SerialJointNode() : Node("serial_joint_node"),count_(0)
  {
    // Set the serial port name
    const char *serialPortName = "/dev/ttyACM0";
    // Set the baud rate
    const LibSerial::BaudRate baudRate = LibSerial::BaudRate::BAUD_115200;

    arduino1.Open(serialPortName);
    arduino1.SetBaudRate(baudRate);

    // test whehter the port opens successfully
    if (!arduino1.IsOpen())
    {
      std::cerr << "Error opening serial port." << std::endl;
    }
    // wait 1 second for Arduino to reboot
    std::this_thread::sleep_for(std::chrono::milliseconds(1000));

    // set QOS
    auto qos = rclcpp::QoS(rclcpp::KeepLast(1), rmw_qos_profile_sensor_data);

    // the publisher will publish the current joint states
    publisher_ = this->create_publisher<communication::msg::JointState>("hardware/joint_states", qos);
    publisher_debug = this->create_publisher<communication::msg::JointState>("hardware/joints_desired_debug", qos);
    // the subscriber will listen to topic trajectory and send the data to arduino
    subscription_ = this->create_subscription<communication::msg::JointState>("hardware/actuators_cmds", qos, std::bind(&SerialJointNode::topic_callback, this, _1));
    // setup the timer and the callback function for sending and reading data from arduino
    timer_ = this->create_wall_timer(std::chrono::milliseconds(5), std::bind(&SerialJointNode::timer_callback, this));

    // initiate the filter 
    // set all the coefficients
    int i; 
    for(i=0;i<4;i++){
      this->filters[i].set_coeff(0.1);
    }
  }

private:
  void topic_callback(communication::msg::JointState msg)
  {
    // store the desired trajectory values
    joint_desired[0] = msg.position[0];
    joint_desired[1] = msg.position[1];
    joint_desired[2] = msg.position[2];
    joint_desired[3] = msg.position[3];
    // store the desired torque compensation
    joint_desired[4] = msg.torque[0];
    joint_desired[5] = msg.torque[1];
    joint_desired[6] = msg.torque[2];
    joint_desired[7] = msg.torque[3];
  }

  void timer_callback()
  { 
    // test sine wave following
    // this->generateSinWaveDebug(); 
    
    this->sendingDataEncode(this->joint_desired, this->joint_desired_encoded);
    this->sendingAndReadingData();
    this->receivedDataDecode(this->receivedData, this->joint_real);

    // read joint real position and filter to get the joint real velocity
    int i; 
    for(i = 0;i<4;i++){
      this->joint_states.position[i] = this->joint_real[i];
      this->joint_states.velocity[i] = this->filters[i].calc(((this->joint_real[i] - this->joint_prev[i])/this->dt));
    }

    // check sanity of the data received
    // theoretically the hip joints should not go beyond 90 degree, knee joints should not go beyond 150 degree 
    
    for(i=0; i<4; i++){
      if(this->joint_real[i] >2.5 || this->joint_real[i] < -2.5){ 
        this->joint_real[i] = this->joint_prev[i]; // use the previous step data 
        this->joint_states.position[i] = this->joint_prev[i]; 
        this->joint_states.velocity[i] = this->joint_velo_prev[i]; 
      }
    }
    // store the positions 
    for(i=0;i<4;i++){
      this->joint_prev[i] = this->joint_real[i];
      // store the velocities to deal with wrong data sending and receiving
      this->joint_velo_prev[i] = this->joint_states.velocity[i]; 
    }

    this->publisher_->publish(this->joint_states);
    this->publisher_debug->publish(this->joint_debug);
  }

  void receivedDataDecode(uint8_t *receivedData, double *joint_states_temp)
  {
    // this method decode the received byte data into double and publish it
    int i;
    uint8_t j1[2];
    uint8_t j2[2];
    uint8_t j3[2];
    uint8_t j4[2];
    for (i = 0; i < 2; i++)
    {
      j1[i] = receivedData[i];
      j2[i] = receivedData[i + 2];
      j3[i] = receivedData[i + 4];
      j4[i] = receivedData[i + 6];
    }
    int16_t temp[8];
    temp[0] = (j1[0] << 8) | j1[1];
    temp[1] = (j2[0] << 8) | j2[1];
    temp[2] = (j3[0] << 8) | j3[1];
    temp[3] = (j4[0] << 8) | j4[1];

    joint_states_temp[0] = temp[0] / 100.0;
    joint_states_temp[1] = temp[1] / 100.0;
    joint_states_temp[2] = temp[2] / 100.0;
    joint_states_temp[3] = temp[3] / 100.0;
  }

  void sendingDataEncode(double *joint_desired, uint8_t *joint_desired_encoded)
  {
    // this method encode the subscribed data from double to char
    int i;
    int16_t temp[8];
    for (i = 0; i < 8; i++)
    {
      temp[i] = (int16_t)(float_round(joint_desired[i]) * 100); // round to two digits, times 100 and type cast to int16_t
    }
    for (i = 0; i < 8; i++)
    {
      joint_desired_encoded[i * 2 + 1] = temp[i] & 0x00ff;  // get the byte in 0x0b
      joint_desired_encoded[i * 2] = (temp[i] >> 8) & 0xff; // get the byte in 0xb0
    }
  }

  float float_round(double var)
  {
    // this function rounds a floating point number to two decimal accuracies
    float value = (int)(var * 100);
    return (float)value / 100;
  }

  void sendingAndReadingData()
  {
    // this method frst read serial port to receive joint states and write serial port to send desired joint states to Arduino
    static bool recInProgress = false;
    uint8_t startMarker = 's';
    uint8_t endMarker = 'e';
    static unsigned char ndx = 0;
    char temp;
    while (this->arduino1.IsDataAvailable())
    {
      arduino1.ReadByte(temp);
      if (recInProgress == true)
      {
        if (temp != endMarker)
        {
          receivedData[ndx] = temp;
          ndx++;

          if (ndx >= 16)
          {
            ndx = 16 - 1;
          }
        }
        else
        {
          recInProgress = false;
          ndx = 0;
        }
      }
      else if (temp == startMarker)
      {
        recInProgress = true;
        ndx = 0;
      }
    }
    int i;
    arduino1.WriteByte('s');
    for (i = 0; i < 16; i++)
    {
      arduino1.WriteByte(this->joint_desired_encoded[i]);
    }
    arduino1.WriteByte('e');
  }

  void generateSinWaveDebug(){
    auto now = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_all_seconds = now - past;
    //std::cout<< elapsed_all_seconds.count() << std::endl; 
    double tem = 0.5*sin(5 * 3.1415 * elapsed_all_seconds.count());
    int i = 0; 
    for(i=0;i<4;i++){
      joint_desired[i] = tem;
      //std::cout << joint_desired[i] <<std::endl; 
      joint_debug.position[i] = tem;
    }
  }
  
  // instantiate the publisher, subscriber and timer
  rclcpp::TimerBase::SharedPtr timer_;
  rclcpp::Publisher<communication::msg::JointState>::SharedPtr publisher_;
  rclcpp::Publisher<communication::msg::JointState>::SharedPtr publisher_debug;
  rclcpp::Subscription<communication::msg::JointState>::SharedPtr subscription_;
  size_t count_;
  std::chrono::time_point<std::chrono::system_clock> past = std::chrono::system_clock::now();

  // build a serial communication object
  SerialPort arduino1;

  uint8_t receivedData[20];
  double joint_real[8] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double joint_prev[4] = {0.0, 0.0, 0.0, 0.0};
  double joint_velo_prev[4] = {0.0,0.0,0.0,0.0};
  double joint_desired[8] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  uint8_t joint_desired_encoded[16] = {0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0};
  double dt = 0.005;

  // instantiate the filter for the joint velocities
  LowPassFilter filters[4];

  // instantiate the message to be published
  communication::msg::JointState joint_states;
  communication::msg::JointState joint_debug;
};
int main(int argc, char *argv[])
{
  rclcpp::init(argc, argv);
  rclcpp::spin(std::make_shared<SerialJointNode>());
  rclcpp::shutdown();
  return 0;
}