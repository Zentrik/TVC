function y = vdpMeasurementNonAdditiveNoiseFcnTranslationOnly(x, v, u)
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
ThrustAcceleration = u(4);
COTCOM = u(5);
q = u(6:9);

Ve = x(1:3);
Xe = x(4:6);

%a_imeas = [0;0;Thrust] + [Mb(1:2);0] / COTCOM + quat2dcm([q(1), q(2), q(3), q(4)]) * [0;0;-9.80655];
a_imeas = [0;0;ThrustAcceleration] + [Mb(1:2);0] / COTCOM;
a_meas = a_imeas + v(1:3);

barometer = Xe(3) + v(4);

y = [a_meas; barometer];
end