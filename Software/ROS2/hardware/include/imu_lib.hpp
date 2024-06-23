/*
    @file: imu_lib.hpp
        This file contains classes for interfacing with IMUs
    
    @Author:
        wjd123ap
        This file is originally open-sourced in ICM20948_ROS_node repo (https://github.com/wjd123ap/ICM20948_ROS_node/blob/main/src/sensor_set/include/sensor_set/IMU_i2c.h)
    
    @Modifier: 
        Yuhao Huang
*/



#include "imu_i2c.h"
#include <iostream>
#include <thread>
uint64_t micro_time(void);
#define GYRO_SCALE_FACTOR 131.0f
#define ACCEL_SCALE_FACTOR 16384.0f
using namespace std;
class ICM20948
{
    public:
        ICM20948();
        void imuDataGet();
        void I2Cinit(int i2c_devnum);
        void MagRead(int16_t& s16X, int16_t& s16Y, int16_t& s16Z);
        void Close(void);
        std::thread imuDataGet_thread(void);
        float halfT;
        float ex_i;
        float ey_i;
        float ez_i;
        float Kp;
        float Ki;
        float q0,q1,q2,q3;
        int *fd_address;
        int sensor_address = 0x68; 
        IMU_QUATERNION_DATA stQuaternion;
        IMU_ST_SENSOR_DATA stGyroRawData;
        IMU_ST_SENSOR_DATA stAccelRawData;
        IMU_ST_SENSOR_DATA stMagRawData;
        IMU_ST_ANGLES_DATA stAngles;
    private:
        
        state_check Check(void);
        void Init(void);
        void GyroRead(int16_t& ps16X, int16_t& ps16Y, int16_t& ps16Z);
        void AccelRead(int16_t& ps16X, int16_t& ps16Y, int16_t& ps16Z);
        void GyroOffset(void);


        void imuAHRSupdate(float gx, float gy, float gz, float ax, float ay, float az, float mx, float my, float mz);
        void ReadSecondary(uint8_t u8I2CAddr, uint8_t u8RegAddr, uint8_t u8Len, uint8_t *pu8data);
        void WriteSecondary(uint8_t u8I2CAddr, uint8_t u8RegAddr, uint8_t u8data);
        void CalAvgValue(uint8_t& Index, int16_t *AvgBuffer, int16_t& InVal, int32_t& OutVal);
        bool MagCheck(void) ;
        
        uint8_t burst_part=6;
        uint8_t burst_all=12;
        IMU_ST_SENSOR_DATA gstGyroOffset;

        ST_AVG_DATA accel_stAvgBuf[3];
        ST_AVG_DATA gyro_stAvgBuf[3];
        ST_AVG_DATA mag_stAvgBuf[3];

        int fd;


};