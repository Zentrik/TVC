#include <Arduino.h>
#include <Wire.h>
#include <Servo.h>
#include <SimpleKalmanFilter.h>
#include <Adafruit_BMP280.h>

Adafruit_BMP280 bmp; // I2C
const int sealevel = 1027; //current sea pressure

Servo servox; //initialise servos for x and y axis (pitch and yaw)  
Servo servoy;

int16_t roll, pitch, yaw; //16bit signed integers, needs to be 16 bit as that is the output of MPU6050

float rollNew, pitchNew, yawNew; //declare variable that will hold current gyro values
float rollOld = 0.0; //initialise variable that hold old gyro values
float pitchOld = 0.0;
float yawOld = 0.0;
double rollAngle, pitchAngle, yawAngle; //declare variable that holds current angle, integrated from angular velocity measurement from gyro

float rollCalibrationValue = -164.63; //offset gyro readings so when no movement reading is near 0
float pitchCalibrationValue = -20.63;
float yawCalibrationValue = 0.0;

SimpleKalmanFilter rollKalmanFilter(1, 1, 0.01); //reduces noise in gyro readings
SimpleKalmanFilter pitchKalmanFilter(1, 1, 0.01);
SimpleKalmanFilter yawKalmanFilter(1, 1, 0.01);

#define LED_PIN PC13 // (Arduino is 13, Teensy is 11, Teensy++ is 6) PC13 is 32
bool blinkState = true; //holds current led status
#define gyro_address 0x68

unsigned long currentMillis; //declare variable that stores the current time
unsigned long previousMillis = 0; //initialise variable that stores time since last PID loop and servo update
unsigned long previousMillisGyro = 0; //initialise variable that stores time since last gyro read
unsigned long previousMillisBaro = 0; //initialise variable that stores time since last barometer read
const int interval = 20;    //20ms = 1/50Hz, frequency of servo
const int intervalGyro = 4; //250Hz, gyro read rate
const int intervalGyroAverage = interval; //gyro reading average rate
const int intervalBaro = 500; //barometer read rate

double MOI = 0.09010865521; //time period squared * string-COM squared * weight/ 4pi squared / string length
int Thrust = 25; //current thrust of motor
int MaxThrust = 25; //maximum thrust of motor

//PID Values, tune with max thrust
const double PIDTVCAngle_P = 0.0444910192906637; 
const double PIDTVCAngle_D = 0.00251899563332339;
const double PIDTVCAngle_I = 0.170401833942346;
const double PIDTVCAngle_N = 91.5676010079333;

//PID Values, maybe they are constants and never need changing I dont know
const int Constant_Value = 1;
const float TSamp_WtEt = 0.01;
const int FilterDifferentiatorTF_InitialS = 0;
const float Integrator_gainval = 0.1;
const int Integrator_IC = 0;
const int FilterDifferentiatorTF_NumCoef0 = 1;
const int FilterDifferentiatorTF_NumCoef1 = -1;

// PID Values for Y axis
/* InitializeConditions for DiscreteTransferFcn: '<S2>/Filter Differentiator TF' */
double FilterDifferentiatorTF_statesy = FilterDifferentiatorTF_InitialS;
/* InitializeConditions for DiscreteIntegrator: '<S1>/Integrator' */
double Integrator_DSTATEy = Integrator_IC;
double rtb_FilterDifferentiatorTFy;
double rtb_Sumy;
double Integrator_tmpy;
double Integratory;
double y = 0; //how much to rotate servo by


// PID Values for Y axis
double FilterDifferentiatorTF_statesx = FilterDifferentiatorTF_InitialS;
double Integrator_DSTATEx = Integrator_IC;
double rtb_FilterDifferentiatorTFx;
double rtb_Sumx;
double Integrator_tmpx;
double Integratorx;
double x = 0; //how much to rotate servo by

//Function declarations
void gyro();
void gyroAverage();
void PIDy(double errorAngle);
void PIDx(double errorAngle);

// ================================================================
// ===                      INITIAL SETUP                       ===
// ================================================================

void setup() {
  Wire.setClock(400000); //Set the i2d clock speed to 400khz
  Wire.begin();          //Start the I2C as master

  Serial.begin(115200); //Start serial connection at 11520bps
  delay(250);           //wait 250ms so gyro can start
  Serial.println(F("Begun"));

  Wire.beginTransmission(gyro_address); //start communication with MPU6050
  Wire.write(0x6B);                     //Start writing at register (PWR_MGMT_1)
  Wire.write(0x00);                     //Set register 0x6B to 0 (wakes up MPU6050)
  Wire.endTransmission();               //terminate connection
  
  Wire.beginTransmission(gyro_address); //Start communication with the MPU-6050.
  Wire.write(0x1B);                     //We want to write to the GYRO_CONFIG register (1B hex).
  Wire.write(0x08);                     //Set the register bits as 00001000 (500dps full scale).
  Wire.endTransmission();               //End the transmission with the gyro.
  
  Wire.beginTransmission(gyro_address); //Start communication with the MPU-6050.
  Wire.write(0x1A);                     //We want to write to the CONFIG register (1A hex).
  Wire.write(0x03);                     //Set the register bits as 00000011 (Set Digital Low Pass Filter to ~43Hz).
  Wire.endTransmission();               //End the transmission with the gyro.

  pinMode(LED_PIN, OUTPUT); //setup LED as output

  //servop.attach(PA0); // attach servos
  //servoy.attach(PA1);

  if (!bmp.begin()) {
    Serial.println(F("Could not find a valid BMP280 sensor, check wiring!"));
    while (1);
  }

  /* Default settings from datasheet. */
  bmp.setSampling(Adafruit_BMP280::MODE_NORMAL,     /* Operating Mode. */
                  Adafruit_BMP280::SAMPLING_X2,     /* Temp. oversampling */
                  Adafruit_BMP280::SAMPLING_X16,    /* Pressure oversampling */
                  Adafruit_BMP280::FILTER_X16,      /* Filtering. */
                  Adafruit_BMP280::STANDBY_MS_500); /* Standby time. */
}

void loop() {
  currentMillis = millis(); //set variable to current time in milliseconds

  if (currentMillis - previousMillisGyro >= intervalGyro) {
    previousMillisGyro = currentMillis; //so currentMillis - previousMillisGyro is now 0
    gyro(); //read gyro values
  }

  if (currentMillis - previousMillis >= interval) {
    previousMillis = currentMillis; // save the last time you blinked the LED
    
    gyroAverage();
    PIDx(rollAngle);        //pass x axis rotation into pid function
    PIDy(pitchAngle);       //pass y axis rotation into pid function

    Serial.print(x);            //pid output for x axis
    Serial.print("\t");
    Serial.print(rollAngle,5);  //x axis rotation
    Serial.print("\t");
    Serial.print(y);            //pid output for y axis
    Serial.print("\t");
    Serial.print(pitchAngle,5); //y axis rotation
    Serial.print("\n");

    //servop.write(x);
    //servoy.write(y);
  }

  if (currentMillis - previousMillisBaro >= intervalBaro) {
    previousMillisBaro = currentMillis;
    //Serial.println(bmp.readAltitude(sealevel)); //read altitude from barometer using current sea pressure

    if (blinkState) { //if true turn off, if false turn on
      digitalWrite(LED_PIN, LOW);
    } else {
      digitalWrite(LED_PIN, HIGH);
    }
    blinkState = !blinkState; //invert blinkstate
  }
}

void gyro() {
  Wire.beginTransmission(gyro_address);   //start communication with MPU6050
  Wire.write(0x43);                       //Send byte 0x43 to indicate startregister
  Wire.endTransmission();                 //terminate connection
  Wire.requestFrom(gyro_address, 6);      //request 6 bytes from MPU6050
  roll = Wire.read() << 8 | Wire.read();  //shift high byte left and add low and high byte to X
  pitch = Wire.read() << 8 | Wire.read(); //shift high byte left and add low and high byte to Y
  yaw = Wire.read() << 8 | Wire.read();   //shift high byte left and add low and high byte to Z

  rollNew += (roll - rollCalibrationValue) / 65.5;      //65.5 = 1 deg/sec (check the datasheet of the MPU-6050 for more information).
  pitchNew += (pitch - pitchCalibrationValue) / 65.5;   //to average them
}

void gyroAverage() {
  rollNew = rollKalmanFilter.updateEstimate(rollNew / (intervalGyroAverage / intervalGyro));    //average gyro readings and then pass into kalman filter
  rollAngle += (rollNew + rollOld) * intervalGyroAverage / 2000;        //trapezoidal integral

  if (rollAngle > 180) {            //constrain angle to be between -180 to + 180
    rollAngle -= 360;               //if angle is 190 it becomes -170
  } else if (rollAngle < -180) {
    rollAngle += 360;               //if angle is -190 it becomes 170
  }
  rollOld = rollNew;
  rollNew = 0;

  pitchNew = pitchKalmanFiter.updateEstimate(pitchNew / (intervalGyroAverage / intervalGyro);   //average gyro readings and then pass into kalman filter
  pitchAngle += (pitchNew + pitchOld) * intervalGyroAverage / 2000;     //trapezoidal integral

  if (pitchAngle > 180) {           //constrain angle to be between -180 to + 180
    pitchAngle -= 360;              //if angle is 190 it becomes -170
  } else if (pitchAngle < -180) {
    pitchAngle += 360;              //if angle is -190 it becomes 170
  }
  pitchOld = pitchNew;
  pitchNew = 0;
}

void PIDy(double errorAngle) {
    errorAngle = errorAngle / 180 * PI; //degrees to rads
    
    rtb_FilterDifferentiatorTFy = PIDTVCAngle_N * TSamp_WtEt;
    rtb_Sumy = 1.0 / (Constant_Value + rtb_FilterDifferentiatorTFy);
    rtb_FilterDifferentiatorTFy = PIDTVCAngle_D * errorAngle - (rtb_FilterDifferentiatorTFy - Constant_Value) * rtb_Sumy * FilterDifferentiatorTF_statesy;
    Integrator_tmpy = Integrator_gainval * (PIDTVCAngle_I * errorAngle);
    Integratory = Integrator_tmpy + Integrator_DSTATEy;
    y = (FilterDifferentiatorTF_NumCoef0 * rtb_FilterDifferentiatorTFy + FilterDifferentiatorTF_NumCoef1 * FilterDifferentiatorTF_statesy) * rtb_Sumy * PIDTVCAngle_N + (PIDTVCAngle_P * errorAngle + Integratory);

    /* Update for DiscreteTransferFcn: '<S2>/Filter Differentiator TF' */
    FilterDifferentiatorTF_statesy = rtb_FilterDifferentiatorTFy; 
    /* Update for DiscreteIntegrator: '<S1>/Integrator' */
    Integrator_DSTATEy = Integrator_tmpy + Integratory;

    if (y < 0.174532925199433 || y > -0.174532925199433) {  // if y > -10 or y < 10. 10 degress = 0.174532925199433 radians. This prevents y > pi/2 when sin(y) would start decreasing
      y = asin(sin(y) * MaxThrust / Thrust);                //scale y so that sin(y) is increased linearly depending on MaxThrust / Thrust. This means the applied moment should stay the same as long as thrust isn't too small that y > 10 or y < -10.
    } 
   
    y = y / PI * 180;       //convert x from radians to degress
    if (y > 10) {           //if greater than 10, x = 10
      y = 10;
    } else if (y < -10) {   //if less than -10, x = -10
      y = -10;
    }  
}

void PIDx(double errorAngle) {
    errorAngle = errorAngle / 180 * PI; //degrees to rads
    
    rtb_FilterDifferentiatorTFx = PIDTVCAngle_N * TSamp_WtEt;
    rtb_Sumx = 1.0 / (Constant_Value + rtb_FilterDifferentiatorTFx);
    rtb_FilterDifferentiatorTFx = PIDTVCAngle_D * errorAngle - (rtb_FilterDifferentiatorTFx - Constant_Value) * rtb_Sumx * FilterDifferentiatorTF_statesx;
    Integrator_tmpx = Integrator_gainval * (PIDTVCAngle_I * errorAngle);
    Integratorx = Integrator_tmpx + Integrator_DSTATEx;
    x = (FilterDifferentiatorTF_NumCoef0 * rtb_FilterDifferentiatorTFx + FilterDifferentiatorTF_NumCoef1 * FilterDifferentiatorTF_statesx) * rtb_Sumx * PIDTVCAngle_N + (PIDTVCAngle_P * errorAngle + Integratorx);

    /* Update for DiscreteTransferFcn: '<S2>/Filter Differentiator TF' */
    FilterDifferentiatorTF_statesx = rtb_FilterDifferentiatorTFx; 
    /* Update for DiscreteIntegrator: '<S1>/Integrator' */
    Integrator_DSTATEx = Integrator_tmpx + Integratorx;

    if (x < 0.174532925199433 || x > -0.174532925199433) { // if x > -10 or x < 10. 10 degress = 0.174532925199433 radians. This prevents x > pi/2 when sin(x) would start decreasing
      x = asin(sin(x) * MaxThrust / Thrust); //scale x so that sin(x) is increased linearly depending on MaxThrust / Thrust. This means the applied moment should stay the same as long as thrust isn't too small that x > 10 or x < -10.
    } 
   
    x = x / PI * 180; //convert x from radians to degress
    if (x > 10) {               //if greater than 10, x = 10
      x = 10;
    } else if (x < -10) {       //if less than -10, x = -10
      x = -10;
    }  
}
