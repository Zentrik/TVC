function y = vdpMeasurementNonAdditiveNoiseFcn(x, v, u)
% vdpMeasurementNonAdditiveNoiseFcn Example measurement function for discrete
% time nonlinear state estimators with non-additive measurement noise.
%
% yk = vdpNonAdditiveMeasurementFcn(xk,vk)
%
% Inputs:
%    xk - x[k], states at time k
%    vk - v[k], measurement noise at time k
%
% Outputs:
%    yk - y[k], measurements at time k
%
% The measurement is the first state with multiplicative noise
%
% See also extendedKalmanFilter, unscentedKalmanFilter

%   Copyright 2016 The MathWorks, Inc.

%#codegen

% The tag %#codegen must be included if you wish to generate code with 
% MATLAB Coder.

Mb = u(1:3);
Acceleration = u(4:6);
I = diag(u(7:9));

w = x(1:3);
q = x(4:7);
Ve = x(8:10);
Xe = x(11:13);

qin = q';
DCM = [qin(:,1).^2 + qin(:,2).^2 - qin(:,3).^2 - qin(:,4).^2 , 2.*(qin(:,2).*qin(:,3) + qin(:,1).*qin(:,4)), 2.*(qin(:,2).*qin(:,4) - qin(:,1).*qin(:,3)) ; ...
    2.*(qin(:,2).*qin(:,3) - qin(:,1).*qin(:,4)), qin(:,1).^2 - qin(:,2).^2 + qin(:,3).^2 - qin(:,4).^2, 2.*(qin(:,3).*qin(:,4) + qin(:,1).*qin(:,2)) ; ...
    2.*(qin(:,2).*qin(:,4) + qin(:,1).*qin(:,3)), 2.*(qin(:,3).*qin(:,4) - qin(:,1).*qin(:,2)), qin(:,1).^2 - qin(:,2).^2 - qin(:,3).^2 + qin(:,4).^2];

w_meas = w + v(1:3) + x(14:16);
%a_imeas = [0;0;Thrust] + [Mb(1:2);0] / COTCOM + quat2dcm([q(1), q(2), q(3), q(4)]) * [0;0;-9.80655];
a_imeas = Acceleration + x(17:19);
if Xe(3) <= 0
    a_imeas = a_imeas + DCM * [0;0;9.80655];
end
a_meas = a_imeas + v(4:6);

barometer = Xe(3) + v(7);

y = [a_meas; w_meas; barometer];
% y = [a_meas; w_meas];
end