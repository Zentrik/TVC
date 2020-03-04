# DEPRECATED PROJECT, WILL BE UPDATED IN 2021 TO WORKING CONTROLLER
# TVC
Quaternion 6DOF.slx is a simulink simulation of the rocket attempting to maintain an upright position during flight. 

Tune2.slx is for tuning PID values.

Main.cpp is an implementation in arduino for the bluepill (stm32f103c8t6).It should work with arduinos, just copy paste the contents into the arduino ide or whichever you prefer

Using [utility header files](https://github.com/adafruit/Adafruit_BNO055/tree/master/utility) from adafuit bno055 library for quaternion class.

# Recommended Watching

[BPS.space](https://www.youtube.com/channel/UCILl8ozWuxnFYXIe2svjHhg)

[Brian Douglas](https://www.youtube.com/user/ControlLectures/videos)

[Understanding PID Control, Part 1: What is PID Control?](https://youtu.be/wkfEZmsQqiA)

[How to read values from mpu6050](https://www.youtube.com/watch?v=ImctYI8hgq4)

[Quaternions by 3Blue1Brown](https://www.youtube.com/watch?v=d4EgbgTm0Bg)

# Recommended Reading
[The Fundamentals of Control Theory](https://www.patreon.com/posts/book-is-now-free-28313078)

[Books Joe Barnard (BPS.space) recommends](https://www.youtube.com/watch?v=BcKL4M5Xod)

[Quaternions](https://folk.uio.no/jeanra/Informatics/QuaternionsAndIMUs.html)

https://folk.ntnu.no/skoge/prost/proceedings/ecc-2013/data/papers/0927.pdf

https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20120014565.pdf

http://ardupilot.org/dev/docs/apmcopter-programming-attitude-control-2.html

http://ardupilot.org/dev/docs/ekf2-estimation-system.html?highlight=quaternion

Instead of trying to estimate the quaternion orientation directly, it estimates an error rotation vector and applies the correction to the quaternion from the inertial navigation equations. This is better when there are large angle errors as it avoids errors associated with linearising quaternions across large angle changes.

    See “Rotation Vector in Attitude Estimation”, Mark E. Pittelkau, Journal of Guidance, Control, and Dynamics, 2003” for details on this approach.
