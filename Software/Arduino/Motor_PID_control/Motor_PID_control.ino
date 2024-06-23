//// This code is about reading the encoder data using interrupts then calculate the PID values and put the PID data into the motors
#include <stdint.h> 
// serial communication variables
int motor_id[4] = {0,1,2,3}; 
const char numChars = 20;
uint8_t receiveddata[numChars];
uint8_t joint_real_encoded[8];
// double joint_states_desired[4] = {0.00,0.00,0.00,0.00};
double joint_states_desired[8] = {0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00}; 
double joint_real[4] = {0.000,0.000,0.000,0.000};

int debug_bit = 0;

// PWM pins
int PWM_pin[4] {7, 8, 11,12};  // 8 7 12 11 
int PWM_pin_cw[4] = {28,31,33,34}; 
int PWM_pin_ccw[4] = {29,30,32,35};
// PWM variables
double actual[4]; 
int direct[4] = {1,1,1,1}; 
double pwm_grav[4] = {0.0,0.0,0.0,0.0};


// Encoder pins
int encoder_pin_cw[4] = {2,18,20,21};
int encoder_pin_ccw[4] = {3,19,14,15}; 
// filter for encoders
class MovingAverageFilter{
public:  
  int window_size; 
  double previous_data[50] = {0}; // the max size of the window is set to 50
  MovingAverageFilter(int size){
    this->window_size = size; 
  } 
  double calc(){
    double output =0; 
    for(int i=0;i<window_size;i++){
      output = output + previous_data[i]; 
    }
    return (output/window_size); 
  }
  void renew(double input){
    for(int i=0;i<window_size;i++){
      previous_data[i]= previous_data[i+1]; 
    }
    previous_data[window_size - 1] = input; 
  }
};
MovingAverageFilter filter[4](2);
// encoder and position readings 
double resolution = 1992.3/(4*2*3.14159);
int count[4] = {0,0,0,0}; 
double position[4] = {0,0,0,0}; 

unsigned long start_time = millis(); 
int count_time_debug = 0; 
class PID{
public:  
  PID(double kp,double kd,double ki):Kp(kp),Kd(kd),Ki(ki){}; 
  double ComputeOutput(){
    // this function will compute the output to the PWM write function
    unsigned long now = millis(); 
    double timeChange = (double)((now - lastTime)); 
    // compute errors
    double error = this->setpoint - this->input; 
    ErrorSum += (error*timeChange); 
    double dErr = (error - lastError) / timeChange; 
    output = Kp*error + Ki*ErrorSum+ Kd*dErr;
    if (output>100){
      output = 100;
    }
    else if (output<-100){
      output = -100; 
    }
    lastError = error; 
    lastTime = now;
    return output;  
  }
  void Setsetpoint(double in){
    this-> setpoint = in;
  }
  void SetRealPosition(double pos){
    this-> input = pos; 
  }
private:   
  double Kp; 
  double Kd; 
  double Ki;  
  double ErrorSum = 0; 
  double lastError; 
  unsigned long lastTime = 0;
  double input = 0; 
  double output = 0;
  double setpoint = 0;   
};
// build up PID controllers
PID pd_controller[4]={PID(100.0,10.0,0),PID(100.0,10.0,0),PID(100.0,10.0,0),PID(100.0,10.0,0)};
void getPosition(){
  int i = 0; 
  for(i=0;i<4;i++){
    joint_real[i] = count[i]/resolution;
    filter[i].renew(joint_real[i]); 
    joint_real[i] = filter[i].calc();
  }  
}
void checkBound(){
  int i = 0;
  for(i=0;i<4;i++){
    if(i == 1 || i == 3){
    if(joint_states_desired[i] > 2.2 ){
      joint_states_desired[i] = 2.2; 
    }
    else if(joint_states_desired[i] < -1.7){
      joint_states_desired[i] = -1.7;
    }
    }
    if(i == 0 || i == 2){
    if(joint_states_desired[i] > 1.7 ){
      joint_states_desired[i] = 1.7; 
    }
    else if(joint_states_desired[i] < -1.7){
      joint_states_desired[i] = -1.7;
    }
    }

  }
}

void setup()
{
  Serial.begin(115200);
  // setup a high PWM frequency
  TCCR1A = 0b10100000;
  TCCR1B = 0b00010001;
  TCCR4A = 0b10100000;
  TCCR4B = 0b00010001;
  ICR1 = 100;
  ICR4 = 100;
  int i;
  for (i=0;i<4;i++){
    // set the encoder pin mode as input 
    pinMode(encoder_pin_cw[i],INPUT_PULLUP);
    pinMode(encoder_pin_ccw[i],INPUT_PULLUP);
    // set the PWM pin mode as output
    pinMode(PWM_pin[i],OUTPUT);
    pinMode(PWM_pin_cw[i],OUTPUT);
    pinMode(PWM_pin_ccw[i],OUTPUT);

  }
  // attach intterupt for encoder pins 
  attachInterrupt(digitalPinToInterrupt(encoder_pin_cw[0]),readEncoder<0>,RISING);
  attachInterrupt(digitalPinToInterrupt(encoder_pin_cw[1]),readEncoder<1>,RISING);
  attachInterrupt(digitalPinToInterrupt(encoder_pin_cw[2]),readEncoder<2>,RISING);
  attachInterrupt(digitalPinToInterrupt(encoder_pin_cw[3]),readEncoder<3>,RISING);
}
template<int j>
void readEncoder(){
/////////////////// method reading two pins /////////////////
  // old_val[j] = new_val[j]; 
  // new_val[j] = digitalRead(encoder_pin_cw[j])*2 + digitalRead(encoder_pin_ccw[j]);
  // out = QEM[old_val[j]*4+new_val[j]]; 
  // count[j] = count[j] + out; 

///////////////// method only read one pin  ///////////////
  int b = digitalRead(encoder_pin_ccw[j]);
  if(b>0){
    count[j]--; // ++ for pull down
  }
  else{
    count[j]++; 
  }
}

void getPosValue_byte(){
    static boolean recvInProgress = false;
    static byte ndx = 0;
    char startMarker = 's'; 
    char endMarker = 'e';
    char rc;
    while (Serial.available() > 0) { 
        rc = Serial.read();
        debug_bit = 1; 
        if (recvInProgress == true) {
            if (rc != endMarker) {
                receiveddata[ndx] = rc;
                ndx++;
                if (ndx >= numChars) {
                    ndx = numChars - 1;
                }
            }
            else {
                recvInProgress = false;
                ndx = 0;
            }
        }         
        if (rc == startMarker) {
            recvInProgress = true;
        }
    }
}

void PWMWrite(int motor_id, double DutyCycle, int direc){
  if (direc == 1){
    digitalWrite(PWM_pin_ccw[motor_id],LOW);
    digitalWrite(PWM_pin_cw[motor_id],HIGH); 
  }
  else if (direc == 0){
    digitalWrite(PWM_pin_cw[motor_id],LOW);
    digitalWrite(PWM_pin_ccw[motor_id],HIGH); 
  }
  double OutPutValue = DutyCycle*255/100;
  analogWrite(PWM_pin[motor_id],int(OutPutValue));
}

void receivedDataDecode(){
  int i;
  uint8_t j1[2];
  uint8_t j2[2];
  uint8_t j3[2];
  uint8_t j4[2];
  uint8_t j5[2];
  uint8_t j6[2];
  uint8_t j7[2];
  uint8_t j8[2];
  
  for (i = 0; i < 2; i++) {
    j1[i] = receiveddata[i];
    j2[i] = receiveddata[i + 2];
    j3[i] = receiveddata[i + 4];
    j4[i] = receiveddata[i + 6];
    j5[i] = receiveddata[i + 8];
    j6[i] = receiveddata[i + 10];
    j7[i] = receiveddata[i + 12];
    j8[i] = receiveddata[i + 14];
  }

  int16_t temp[8];
  temp[0] = (j1[0] << 8) | j1[1];
  temp[1] = (j2[0] << 8) | j2[1];
  temp[2] = (j3[0] << 8) | j3[1];
  temp[3] = (j4[0] << 8) | j4[1];
  temp[4] = (j5[0] << 8) | j5[1];
  temp[5] = (j6[0] << 8) | j6[1];
  temp[6] = (j7[0] << 8) | j7[1];
  temp[7] = (j8[0] << 8) | j8[1];

  joint_states_desired[0] = temp[0] / 100.0;
  joint_states_desired[1] = temp[1] / 100.0;
  joint_states_desired[2] = temp[2] / 100.0;
  joint_states_desired[3] = temp[3] / 100.0;
  joint_states_desired[4] = temp[4] / 100.0;
  joint_states_desired[5] = temp[5] / 100.0;
  joint_states_desired[6] = temp[6] / 100.0;
  joint_states_desired[7] = temp[7] / 100.0; 
}

void torquePWMMap(){
  // this method maps torque compensation values to PWM values using: v_pwm = R/(KT*N)*tau
  // motor params: N = 71.2, R = 1.3043, KT = 0.014
  double temp[4] = {0.0,0.0,0.0,0.0};
  int i; 
  for(i=0;i<4;i++){
    temp[i] = 1.3043/(0.014*71.2)*joint_states_desired[i+4];
    pwm_grav[i] = map(temp[i],-12.0,12.0,-100.0,100.0); 
  }
}

void sendingDataEncode() {
  // this method encode the subscribed data from double to char
  int i;
  int16_t temp[4];
  for (i = 0; i < 4; i++) {
    temp[i] = (int16_t)(float_round(joint_real[i]) * 100);  // round to two digits, times 100 and type cast to int8_t
  }
  for (i = 0; i < 4; i++) {
    joint_real_encoded[i * 2 + 1] = temp[i] & 0x00ff;   // get the byte in 0x0b
    joint_real_encoded[i * 2] = (temp[i] >> 8) & 0xff;  // get the byte in 0xb0
                                                           // cout<<hex<<(unsigned char)joint_desired_encoded[i*2]<<endl;
  }
}

float float_round(double var) {
  // this function rounds a floating point number to two decimal accuracies
  float value = (int)(var * 100);
  return (float)value / 100;
}

long int prev_val;
double generate_sine_wave_debug(){
  // step input
  // return 3.14; // rotate 180 degree
  
  // sine wave
  double period = 0.01;
  unsigned long t = millis();  
  return 3*sin((t-start_time)/period*6.28318);

  // square wave
  // if(t-start_time > 1000){
  //   if(prev_val == 0){
  //     prev_val = 3;
  //     return 3;
  //   }
  //   else{
  //     prev_val = 0;
  //     return 0;
  //   }
  //   start_time = t; 
  // }

}

void loop(){  
  // generate a sine wave for testing
  //double trajectory_debug = generate_sine_wave_debug(); 

  // read desired trajectory from Raspberry Pi 
  getPosValue_byte(); 
  receivedDataDecode();
  debug_bit == 0;
  getPosition(); 
  // send the real joint values
  sendingDataEncode();
  int i;
  Serial.write('s');
  for(i=0;i<8;i++){
    Serial.write(joint_real_encoded[i]);
  }
  Serial.write('e');

  // do PD control for motors 
  checkBound();
  for(i=0;i<4;i++){
    pd_controller[i].Setsetpoint(joint_states_desired[i]); 
    pd_controller[i].SetRealPosition(joint_real[i]);
    torquePWMMap();
    double PWM = pd_controller[i].ComputeOutput()+pwm_grav[i]; 
    if(PWM < 0){
      PWM = -PWM; 
      direct[i] = 0; 
    }
    else{
      direct[i] = 1; 
    }
    PWMWrite(i,PWM,direct[i]); 
 }
  long int b = micros();
}
