#include <iostream>
#include <wiringPi.h>

#include "rclcpp/rclcpp.hpp"
//#include "planar_biped_message_package/msg/contact_sensor.hpp"
#include "communication/msg/contact_state.hpp"

#define HIGH 1
#define LOW 0 

using namespace std::chrono_literals;

/*************************This node is used to read contact sensor data****************************/
enum class contact_state
{
    Waiting,
    Detected,
    Pressed,
    Released
};

class ContactPublisher : public rclcpp::Node
{
public:
    ContactPublisher() : Node("contact_info_publisher"), count_(0)
    {
        ////////// WiringPi config
        wiringPiSetup();
        // // set the contact sensor GPIO pins to 20 and 21
        pinMode(25, INPUT);
        pinMode(27, INPUT);
        pullUpDnControl(25,PUD_DOWN); //GPIO 16
        pullUpDnControl(27,PUD_DOWN); //GPIO 26

        auto qos = rclcpp::QoS(rclcpp::KeepLast(1), rmw_qos_profile_sensor_data);

        // initialize the publisher with a message type and topic name, queue size
        publisher_ = this->create_publisher<communication::msg::ContactState>("hardware/ContactState",qos);
        // setup the timer and the callback function 500ms measn one second run the callback function two times
        timer_ = this->create_wall_timer(5ms, std::bind(&ContactPublisher::timer_callback, this));
    }

private:
    void timer_callback()
    {
        // first read whether the button is activated
        // wiringPi
        state_temp[0] = digitalRead(contact_pin1); 
        state_temp[1] = digitalRead(contact_pin2); 
        int i;
        for (i = 0; i < 2; i++)
        {
            switch (states[i])
            {
            case (contact_state::Waiting):
                //RCLCPP_INFO(this->get_logger(),"waiting");
                if (state_temp[i] == HIGH)
                {
                    states[i] = contact_state::Detected;
                    // RCLCPP_INFO(this->get_logger(),"shift to detected");
                    break; 
                }
                // add contact data to be published
                if(i == 0){
                    this->contact_data.right_contact.data = false; 
                }
                else if(i == 1){
                    this->contact_data.left_contact.data = false; 
                }
                break;
            case (contact_state::Detected):
                //RCLCPP_INFO(this->get_logger(),"detected");
                if (state_temp[i] == HIGH)
                {
                    state_count_up[i]++;
                }
                if (state_count_up[i] > 3)
                {
                    states[i] = contact_state::Pressed;
                    // RCLCPP_INFO(this->get_logger(),"shift to pressed");
                    break;
                }
                if (state_temp[i] == LOW){
                    states[i] = contact_state::Waiting; 
                    break; 
                }
                // add contact data to be published
                if(i == 0){
                    this->contact_data.right_contact.data = false; 
                }
                else if(i == 1){
                    this->contact_data.left_contact.data = false; 
                }
                break;
            case (contact_state::Pressed):
                //RCLCPP_INFO(this->get_logger(),"pressed");
                if (state_temp[i] == LOW){
                    states[i] = contact_state::Released; 
                    break;
                }
                // add contact data to be published
                if(i == 0){
                    this->contact_data.right_contact.data = true; 
                }
                else if(i == 1){
                    this->contact_data.left_contact.data = true; 
                }
                break;
            case (contact_state::Released): 
                //RCLCPP_INFO(this->get_logger(),"released");
                if(state_temp[i] == LOW){
                    state_count_down[i]++; 
                }
                if (state_count_down[i] > 3)
                {
                    states[i] = contact_state::Waiting;
                    break;
                }
                if (state_temp[i] == HIGH)
                {
                    states[i] = contact_state::Pressed;
                    break;
                }
                // add contact data to be published
                if(i == 0){
                    this->contact_data.right_contact.data = false; 
                }
                else if(i == 1){
                    this->contact_data.left_contact.data = false; 
                }
                break;
            }
        }
        rclcpp::Time now = this->get_clock()->now(); 
        contact_data.header.stamp = now; 

        this->publisher_->publish(contact_data);
    }
    rclcpp::TimerBase::SharedPtr timer_;
    rclcpp::Publisher<communication::msg::ContactState>::SharedPtr publisher_;
    size_t count_;
    communication::msg::ContactState contact_data;
    // the contact counter
    contact_state states[2] = {contact_state::Waiting, contact_state::Waiting};
    int state_count_up[2] = {0, 0};
    int state_count_down[2] = {0, 0};
    int state_temp[2] = {0, 0};
    // GPIO pins
    // wiringPi 
    int contact_pin1 = 27;
    int contact_pin2 = 25; 
};
int main(int argc, char *argv[])
{
    rclcpp::init(argc, argv);
    rclcpp::spin(std::make_shared<ContactPublisher>());
    rclcpp::shutdown();
    return 0;
}