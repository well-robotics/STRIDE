#include "DShot.h"

uint16_t bound_target(int target){
    if(target < 48){
        target = 48; 
    }
    else if(target > 2047){
        target = 2047; 
    }
    return target; 
}

uint16_t get_throttle(int target, int throttle_now){
  if(throttle_now < 48){
    throttle_now = 48;
  }
  else if(throttle_now > target){
    throttle_now --; 
  }
  else if(throttle_now < target){
    throttle_now ++; 
  }
  return throttle_now;
}