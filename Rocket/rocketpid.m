Kp = 1;
Ki = 0.1;  
Kd = 3;
Tf = 0.5;
C = pid(Kp,Ki,Kd,Tf);

angularAccelerationX = 25 * C * 0.4/0.09; 
thetaX = cumtrapz(cumtrapz(angularAccelerationX));
