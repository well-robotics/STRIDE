#include "DShot.h"
#include "utils.h"
#include <XBOXUSB.h> 
#ifdef dobogusinclude
#include <spi4teensy3.h>
#endif
#include <SPI.h>

// USB object
USB usb;
XBOXUSB Xbox(&usb);

int16_t speed_des_right = 0;
int16_t speed_des_up = 0;
int LT_state = 0; // if it is pressed, the upper quadcopter will push
int RT_state = 0; // if it is pressed, the front quadcopter will push
int reset = 0; // used to reset the whole system
double force_side = 300.0; // unit in g
double force_up = 500.0; // unit in g
/*

redefine DSHOT_PORT if you want to change the default PORT

Defaults
UNO: PORTD, available pins 0-7 (D0-D7)
Leonardo: PORTB, available pins 4-7 (D8-D11)

e.g.
#define DSHOT_PORT PORTD
*/
DShot esc1(DShot::Mode::DSHOT300);
DShot esc2(DShot::Mode::DSHOT300);

uint16_t throttle1 = 0;
uint16_t target1 = 0;
uint16_t throttle2 = 0;
uint16_t target2 = 0;

int j = 0;

void setup() {
  Serial.begin(115200);
  Serial.print(F("\r\nopen"));
  // check the USB
  if (usb.Init() == -1) {
    Serial.print(F("\r\nOSC did not start"));
    while (1); //if do not start, then halt
  }
  usb.Task();
  while (!Xbox.Xbox360Connected) {
    delay(1);
    usb.Task();
    Serial.print(F("\r\nwait for Xbox"));
  }

  // Use pin 5 and 6
  esc1.attach(5);  
  esc2.attach(6); 
  // unlock the motor controller
  esc1.setThrottle(throttle1);
  esc2.setThrottle(throttle2);
  delay(3000);
}

void loop() {
  /////////// 1. get desired speed from joystick
  usb.Task();
  if (Xbox.Xbox360Connected) {
    speed_des_right = Xbox.getAnalogHat(LeftHatY);
    speed_des_up = Xbox.getAnalogHat(RightHatY);
    if(LT_state == 0){
      LT_state = Xbox.getButtonPress(LT);
    }
    if(RT_state == 0){
      RT_state = Xbox.getButtonPress(RT);
    }
    if(Xbox.getButtonClick(LB)){
      speed_des_right = 0;
      speed_des_up = 0;
      LT_state = (int)0; // if it is pressed, the upper quadcopter will push
      RT_state = (int)0; // if it is pressed, the front quadcopter will push
    }
  } 
  if(speed_des_right < 0){
    speed_des_right = abs(speed_des_right + 1);
  }
  if(speed_des_up < 0){
    speed_des_up = abs(speed_des_up + 1);
  }
  uint16_t target_right = map(speed_des_right,0,32767,48,2047);
  uint16_t target_up = map(speed_des_up,0,32767,48,2047);
  //Serial.println(target);
  ////////// 2. drive the motors to the desired speed 
  // throttle should be in the range below:
  // half: 0x10000000000  
  // full: 0x11111111111
  target1 = target_right;
  target2 = target_up; 
  // get the desired throttle and send
  target1 = bound_target(target1);
  target2 = bound_target(target2); 
  throttle1 = get_throttle(target1, throttle1);
  throttle2 = get_throttle(target2, throttle2); 
  if( LT_state >0){ 
    throttle1 = map(force_side,0,1271.51,48,2047); 
  }
  if( RT_state >0){
    throttle2 = map(force_up,0,1271.51,48,2047);
  }
  esc1.setThrottle(throttle1);
  esc2.setThrottle(throttle2);
  // printout the estimated thrust force
  if (j == 100){
    double thrust_right = (target_right - 48.0)/(2047.0-48.0)*1271.51; 
    double thrust_up = (target_up - 48.0)/(2047.0-48.0)*1271.51;
    Serial.println("Thrust force right"); 
    // Serial.println(LT_state);
    Serial.println(thrust_right); 
    Serial.println("Thrust force up"); 
    Serial.println(thrust_up);
    j = 0;
  }
  j++;
  delay(1);
}
