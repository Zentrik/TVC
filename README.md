Rocket/roll_control.slxc should be the simulator and a controller and state estimator can be generated using an example from Rocket/StateSpaceController.m.

The gain scheduled controllers are incorrect as I linearised about points which are not in equilbrium but did not account for the affine terms.
If I get the time, I will add an iLQR, try to implement partial feedback linearisation and fix the guidance MPC.

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

[My paper](Paper/Paper.pdf). The gain scheduling section is wrong as I didn't linearise about equilibrium points and didn't add an affine term to the state space model to account for this. Also, the appendix is not complete. However, hopefully it is helpful and the references are also recommended to be read.

[The Fundamentals of Control Theory](https://www.patreon.com/posts/book-is-now-free-28313078)

[Books Joe Barnard (BPS.space) recommends](https://www.youtube.com/watch?v=BcKL4M5Xod)

[Quaternions](https://folk.uio.no/jeanra/Informatics/QuaternionsAndIMUs.html)
