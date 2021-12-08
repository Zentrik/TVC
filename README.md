[Roll control](Rocket/roll_control.slx) is the full simulator and a state estimator can be generated using an example from [StateSpaceController](Rocket/StateSpaceController.m). A controller can be generate by modifying [affine LQR](Rocket/affine_controller.m) with an example from [State Space Controller](Rocket/StateSpaceController.m).

[State Space Controller](Rocket/StateSpaceController.m) will be cleaned up soon, I don't think the dynamics are right in any of the examples.

Guidance for pinpoint landing is [here](https://github.com/Zentrik/scp_new_problem).

If I get the time, I will try to implement partial feedback linearisation.

<a href="https://www.paypal.com/paypalme/zentriktvc/5/"><img src="https://raw.githubusercontent.com/andreostrovsky/donate-with-paypal/master/blue.svg" height="40"></a>  
If you found my work useful, consider donating. Thanks.

# Recommended Watching

[BPS.space](https://www.youtube.com/channel/UCILl8ozWuxnFYXIe2svjHhg)

[Brian Douglas](https://www.youtube.com/user/ControlLectures/videos)

[Understanding PID Control, Part 1: What is PID Control?](https://youtu.be/wkfEZmsQqiA)

[Control Bootcamp](https://www.youtube.com/playlist?list=PLMrJAkhIeNNR20Mz-VpzgfQs5zrYi085m)

Other videos by brian, brunton, 3b1b etc.

[How to read values from mpu6050](https://www.youtube.com/watch?v=ImctYI8hgq4)

[Quaternions by 3Blue1Brown](https://www.youtube.com/watch?v=d4EgbgTm0Bg)

# Recommended Reading

[My paper](Paper/Paper.pdf). The gain scheduling section is wrong as I didn't linearise about equilibrium points and didn't add an affine term to the state space model to account for this. Also, the appendix is not complete. However, hopefully it is helpful and the references are also recommended to be read.

[The Fundamentals of Control Theory](https://www.patreon.com/posts/book-is-now-free-28313078)

[Books Joe Barnard (BPS.space) recommends](https://www.youtube.com/watch?v=BcKL4M5Xod)

[Quaternions](https://folk.uio.no/jeanra/Informatics/QuaternionsAndIMUs.html), this link seems to be dead.

[Satellite Dynamics and Control in a Quaternion
Formulation](https://orbit.dtu.dk/files/98594729/Satdyn_mb_2010f.pdf)

[Great introduction to control theory](controls-in-frc.link)

## Dynamics

The kinematics of a 6 DOF rigid body are described [here](https://mathworks.com/help/aeroblks/6dofeulerangles.html#mw_2f302a65-767b-4836-81d3-8d9423421b84) and the quaternion derivative is [here](https://mathworks.com/help/aeroblks/customvariablemass6dofquaternion.html) (this is just standard 1/2 * quaternion * angular velocity (in body coordinates)).

Stevens, Brian, and Frank Lewis. Aircraft Control and Simulation, 2nd ed. Hoboken, NJ: John Wiley & Sons, 2003, explains the equations used in the matlab 6dof block nicely.

I think it's quite important to understand reference frames and coordinate systems so reading the book above or the other one [here](https://mathworks.com/help/aeroblks/6dofeulerangles.html#References) should be helpful in that endeavour.

### My notes on 6dof block

The block takes in the forces, moments and Inertia tensor in the body coordinate system.
V_b is the velocity of the centre of mass w.r.t. (I think the flat earth reference frame) expressed in the body coordinate system whilst V_e is expressed in the flat earth coordinate system.
A_bb is the derivative of V_b wrt to the body reference frame whilst A_be is wrt to the flat earth reference frame. 

### Aerodynamics

[OpenRocket](https://github.com/openrocket/openrocket/releases/download/Development_of_an_Open_Source_model_rocket_simulation-thesis-v20090520/Development_of_an_Open_Source_model_rocket_simulation-thesis-v20090520.pdf) describes the aerodynamic forces acting on a model rocket. 
I don't know if I would recommend trying to implement the equation from here, I struggled quite a lot doing that but using debugging tools on the OpenRocket code helped.
It would probably be best to look at other implementations of the barrowman equations/ find the deleted code in this repo where I had it implemented for a rocket with a conical nose cone and no fins.

## Optimal Landing

[SCP Toolbox and links to papers on SCvx, PTR and GuSTO](https://www.malyuta.name/optimization/tooling/2021/07/15/scp-tutorial.html)

[Post on GFOLD](https://tealquaternion.netlify.app/post/gfold-2007/)