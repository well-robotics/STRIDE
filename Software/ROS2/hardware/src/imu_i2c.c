/*
    @file: imu_i2c.c
        This file contains functions for interfacing with i2c devices
    
    @Author:
        wjd123ap
        This file is originally open-sourced in ICM20948_ROS_node repo (https://github.com/wjd123ap/ICM20948_ROS_node/blob/main/src/sensor_set/include/sensor_set/IMU_i2c.h)
    
    @Modifier: 
        Yuhao Huang
*/



#include "imu_i2c.h"

#ifdef __cplusplus
extern "C" {
#endif


int i2cInit(int i2c_devnum)
{

  int fd;
  char i2c_dev[30];
  sprintf(i2c_dev,"/dev/i2c-%d",i2c_devnum);

  if ((fd = open(i2c_dev, O_RDWR)) < 0)
  {
    
    printf("Failed to open the i2c bus.\n");
  }
  
  return fd;
}

void i2cClose(int *fd_address)
{
  close(*fd_address);
}
uint8_t I2C_ReadOneByte(uint8_t DevAddr, uint8_t RegAddr, int *fd_address)
{
  uint8_t u8Ret;
  if (ioctl(*fd_address, I2C_SLAVE, DevAddr) < 0)
  {
    printf("Failed to acquire bus access and/or talk to slave.\n");
    return 0;
  }
  write(*fd_address, &RegAddr,1);
  read(*fd_address, &u8Ret, 1);
  return u8Ret;
}

void I2C_ReadBurstByte(uint8_t DevAddr, uint8_t RegAddr, uint8_t* buffer,uint8_t count, int *fd_address){
  if (ioctl(*fd_address, I2C_SLAVE, DevAddr) < 0)
  {
    printf("Failed to acquire bus access and/or talk to slave.\n");
    return ;
  }
  write(*fd_address, &RegAddr,1);
  read(*fd_address, buffer, count);
  return;
}

void I2C_WriteBurstByte(uint8_t DevAddr, uint8_t RegAddr, uint8_t* value,uint8_t count, int *fd_address){
  int8_t *buf;

  if (ioctl(*fd_address, I2C_SLAVE, DevAddr) < 0)
  {
    printf("Failed to acquire bus access and/or tCan't use SMBus Quick Write command, will skip some addressesalk to slave.\n");
    return;
  }
  buf = malloc(count);
  buf[0] = RegAddr;
  buf[1] = value[0];
  write(*fd_address, buf, count+1);
  free(buf);
  return;
}

void I2C_WriteOneByte(uint8_t DevAddr, uint8_t RegAddr, uint8_t value, int *fd_address)
{
  int8_t *buf;

  if (ioctl(*fd_address, I2C_SLAVE, DevAddr) < 0)
  {
    printf("Failed to acquire bus access and/or tCan't use SMBus Quick Write command, will skip some addressesalk to slave.\n");
    return;
  }
  buf = malloc(2);
  buf[0] = RegAddr;
  buf[1] = value;
  write(*fd_address, buf, 2);
  free(buf);
  return;
}
float invSqrt(float x) 
{
  float halfx = 0.5f * x;
  float y = x;
  
  long i = *(long*)&y;                //get bits for floating value
  i = 0x5f3759df - (i >> 1);          //gives initial guss you
  y = *(float*)&i;                    //convert bits back to float
  y = y * (1.5f - (halfx * y * y));   //newtop step, repeating increases accuracy
  
  return y;
}

#ifdef __cplusplus
}
#endif