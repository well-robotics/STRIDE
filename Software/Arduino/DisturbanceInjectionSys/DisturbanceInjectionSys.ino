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

void setup() {
  Serial.begin(115200);
  
  // check the USB
  if (usb.Init() == -1) {
    Serial.print(F("\r\nOSC did not start"));
    while (1); //if do not start, then halt
  }
  usb.Task();
  while (!Xbox.Xbox360Connected) {
    delay(1);
    usb.Task();
  }

  // Use pin 7 and 11
  esc1.attach(6);  
  esc2.attach(5); 
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
  } 
  if(speed_des_right < 0){
    speed_des_right = abs(speed_des_right + 1);
  }
  if(speed_des_up < 0){
    speed_des_up = abs(speed_des_up + 1);
  }
  uint16_t target_right = map(speed_des_right,0,32767,48,2047);
  uint16_t target_up = map(speed_des_up,0,32767,48,2047);

  target1 = target_right;
  target2 = target_up; 
  // get the desired throttle and send
  target1 = bound_target(target1);
  target2 = bound_target(target2); 
  throttle1 = get_throttle(target1, throttle1);
  throttle2 = get_throttle(target2, throttle2); 
  esc1.setThrottle(throttle1);
  esc2.setThrottle(throttle2);
  delay(1);
}
