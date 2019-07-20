#include <Arduino.h>
#include <Wire.h>
#include <Servo.h>
#include <SimpleKalmanFilter.h>
#include <Adafruit_BMP280.h>

Adafruit_BMP280 bmp; // I2C
const int sealevel = 1027;

Servo servop;
Servo servoy;

int16_t roll, pitch, yaw; //16bit signed integers, needs to be 16 bit as that is the output of MPU6050

float rollNew, pitchNew, yawNew;
float rollOld = 0.0;
float pitchOld = 0.0;
float yawOld = 0.0;
double rollAngle, pitchAngle, yawAngle;

float rollCalibrationValue = -160.27;
float pitchCalibrationValue = -6.71;
float yawCalibrationValue = 60.02;

SimpleKalmanFilter rollKalmanFilter(1, 1, 0.01);
SimpleKalmanFilter pitchKalmanFilter(1, 1, 0.01);
SimpleKalmanFilter yawKalmanFilter(1, 1, 0.01);

#define LED_PIN 32 // (Arduino is 13, Teensy is 11, Teensy++ is 6) PC13 is 32
bool blinkState = true;
#define gyro_address 0x68

unsigned long currentMillis;
unsigned long previousMillis = 0;
unsigned long previousMillisGyro = 0;
unsigned long previousMillisBaro = 0;
const int interval = 20;    //20ms = 1/50Hz
const int intervalGyro = 4; //20ms = 1/50Hz
const int intervalBaro = 500;

double MOI = 0.09010865521; //time period squared * string-COM squared * weight/ 4pi squared / string length
int Thrust = 25;
int MaxThrust = 25;

double Integrator_DSTATE = 0;
double Filter_DSTATE = 0;
double y = 0;
double FilterCoefficient = 0;

double Integrator_DSTATEx = 0;
double Filter_DSTATEx = 0;
double x = 0;
double FilterCoefficientx = 0;

const double Kp = 0.382662534441352;
const double Kd = 0.139903294226043;
const double Ki = 0.206393482611663;
const double N = 36.46824066771;

void gyro();
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

  //servop.attach(PA0);
  //servoy.attach(PA1);

  /*if (!bmp.begin()) {
    Serial.println(F("Could not find a valid BMP280 sensor, check wiring!"));
    while (1);
  }*/

  /* Default settings from datasheet. */
  bmp.setSampling(Adafruit_BMP280::MODE_NORMAL,     /* Operating Mode. */
                  Adafruit_BMP280::SAMPLING_X2,     /* Temp. oversampling */
                  Adafruit_BMP280::SAMPLING_X16,    /* Pressure oversampling */
                  Adafruit_BMP280::FILTER_X16,      /* Filtering. */
                  Adafruit_BMP280::STANDBY_MS_500); /* Standby time. */
}

void loop() {
  currentMillis = millis();

  if (currentMillis - previousMillisGyro >= intervalGyro) {
    previousMillisGyro = currentMillis;
    gyro();
  }

  if (currentMillis - previousMillis >= interval) {
    // save the last time you blinked the LED
    previousMillis = currentMillis;

    /* Serial.print(rollAngle,5);
    Serial.print("\t");
    Serial.print(pitchAngle, 5);
    Serial.print("\t");
    Serial.print(yawAngle, 5); //print gyro values
    Serial.print("\n"); */

    PIDx(rollAngle);
    PIDy(pitchAngle);


    //Serial.print(x);
    //Serial.print("\t");
    Serial.print(rollAngle / PI * 180,5);
    Serial.print("\t");
    //Serial.print(y);
    //Serial.print("\t");
    Serial.print(pitchAngle / PI * 180,5);
    Serial.print("\n");

    //servop.write(x);
    //servoy.write(y);

    if (blinkState) { //if true turn off, if false turn on
      digitalWrite(LED_PIN, LOW);
    } else {
      digitalWrite(LED_PIN, HIGH);
    }
    blinkState = !blinkState; //invert blinkstate
  }

  if (currentMillis - previousMillisBaro >= intervalBaro) {
    previousMillisBaro = currentMillis;
    //Serial.println(bmp.readAltitude(sealevel));
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

  //65.5 = 1 deg/sec (check the datasheet of the MPU-6050 for more information).
  //rollNew = (rollOld * 0.7) + ((((float)roll - rollCalibrationValue)/ 65.5) * 0.3);   //Gyro pid input is deg/sec.
  rollNew = rollKalmanFilter.updateEstimate((roll - rollCalibrationValue) / 65.5);
  pitchNew = pitchKalmanFilter.updateEstimate((pitch - pitchCalibrationValue) / 65.5); //Gyro pid input is deg/sec.
  yawNew = yawKalmanFilter.updateEstimate((yaw - yawCalibrationValue) / 65.5);         //Gyro pid input is deg/sec.

  /* if (abs(rollNew) <= 0.2) {
    rollNew = 0;
  }
  if (abs(pitchNew) <= 0.2) {
    pitchNew = 0;
  }
  if (abs(yawNew) <= 0.2) {
    yawNew = 0;
  } */

  /* Gyro angle calculations
  //T = 0.004 = (1/250Hz), trapezium area is (new + old)/2 * T
  rollAngle += (rollNew + rollOld) * 0.002;    //Calculate the traveled roll angle and add this to the angle_roll variable.
  pitchAngle += (pitchNew + pitchOld) * 0.002; //Calculate the traveled pitch angle and add this to the angle_pitch variable.
  yawAngle += (yawNew + yawOld) * 0.002;*/

  /*0.000001066 = 0.0000611 * (3.142(PI) / 180degr) The Arduino sin function is in radians and not degrees.
  angle_pitch -= angle_roll * sin((float)gyro_yaw * 0.000001066);                  //If the IMU has yawed transfer the roll angle to the pitch angel.
  angle_roll += angle_pitch * sin((float)gyro_yaw * 0.000001066);*/
  //If the IMU has yawed transfer the pitch angle to the roll angel.

  pitchAngle += (pitchNew + pitchOld) * (intervalGyro * PI / 360000);                                    //Calculate the traveled pitch angle and add this to the angle_pitch variable.
  rollAngle += (rollNew + rollOld) * (intervalGyro * PI / 360000);                                 //Calculate the traveled roll angle and add this to the angle_roll variable.
  yawAngle += (yawNew + yawOld) * (intervalGyro * PI / 360000);

  rollOld = rollNew;
  pitchOld = pitchNew;
  yawOld = yawNew;
}

void PIDy(double errorAngle) {
    FilterCoefficient = (Kd * errorAngle - Filter_DSTATE) * N;
    y = (Kp * errorAngle + Integrator_DSTATE) + FilterCoefficient;
    Integrator_DSTATE += Ki * errorAngle * interval / 1000;
    Filter_DSTATE += interval / 1000 * FilterCoefficient;

    y = asin(sin(y) * MaxThrust / Thrust) / PI * 180;
}

void PIDx(double errorAngle) {
    FilterCoefficientx = (Kd * errorAngle - Filter_DSTATEx) * N;
    x = (Kp * errorAngle + Integrator_DSTATEx) + FilterCoefficientx;
    Integrator_DSTATEx += Ki * errorAngle * interval / 1000;
    Filter_DSTATEx += interval / 1000 * FilterCoefficientx;

    x = asin(sin(x) * MaxThrust / Thrust) / PI * 180;
}