# DEPRECATED PROJECT, WILL BE UPDATED IN 2021 TO WORKING GNC SCHEME

Rocket/roll_control.slxc should be the simulator and a controller and state estimator can be generated using an example from Rocket/StateSpaceController.m.

<a href="https://www.paypal.com/paypalme/zentriktvc/5/"><img src="https://raw.githubusercontent.com/andreostrovsky/donate-with-paypal/master/blue.svg" height="40"></a>  
If you found my work useful, conside donating. Thanks.

# Recommended Watching

[BPS.space](https://www.youtube.com/channel/UCILl8ozWuxnFYXIe2svjHhg)

[Brian Douglas](https://www.youtube.com/user/ControlLectures/videos)

[Understanding PID Control, Part 1: What is PID Control?](https://youtu.be/wkfEZmsQqiA)

Other videos by brian, brunton, 3b1b etc.

[How to read values from mpu6050](https://www.youtube.com/watch?v=ImctYI8hgq4)

[Quaternions by 3Blue1Brown](https://www.youtube.com/watch?v=d4EgbgTm0Bg)

# Recommended Reading

My paper, it's not great but its something.

[The Fundamentals of Control Theory](https://www.patreon.com/posts/book-is-now-free-28313078)

[Books Joe Barnard (BPS.space) recommends](https://www.youtube.com/watch?v=BcKL4M5Xod)

[Quaternions](https://folk.uio.no/jeanra/Informatics/QuaternionsAndIMUs.html)

https://folk.ntnu.no/skoge/prost/proceedings/ecc-2013/data/papers/0927.pdf

https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20120014565.pdf

http://ardupilot.org/dev/docs/apmcopter-programming-attitude-control-2.html

http://ardupilot.org/dev/docs/ekf2-estimation-system.html?highlight=quaternion

Instead of trying to estimate the quaternion orientation directly, it estimates an error rotation vector and applies the correction to the quaternion from the inertial navigation equations. This is better when there are large angle errors as it avoids errors associated with linearising quaternions across large angle changes.

    See “Rotation Vector in Attitude Estimation”, Mark E. Pittelkau, Journal of Guidance, Control, and Dynamics, 2003” for details on this approach.
