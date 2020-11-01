ThrustTable = xlsread('F15_Thrust','Sheet1');
mass = 1.0565;
g = 9.80655;
ThrustAcceleration = griddedInterpolant([-1; ThrustTable(:,1); 100], [0; ThrustTable(:,2); 0] / mass);

tstart = t0;
tfinal = 1000;

y0 = [0; 0];
options = odeset('Events',@events,'Refine',4, 'RelTol', 1e-4, 'AbsTol', 1e-8);

yeout = [];

while tstart < 3.45
%     [t,y,~,ye,~] = ode45(@(t,y) ydot(t, y, ThrustAcceleration), [tstart tfinal], y0, options);
    [t,y] = ode45(@(t,y) ~(y(1) <= 0 && ThrustAcceleration(t) <= g) * [y(2); ThrustAcceleration(t) - g], [tstart tfinal], y0, options)

    yeout = [yeout; y(end)];

    y0 = [0;0];
    tstart = t(end);
end
max(abs(yeout))

% function dydt = ydot(t, y, ThrustAcceleration)
%     dydt = ~(y(1) <= 0 && (ThrustAcceleration(t) - 9.80655) <= 0) * [y(2); (ThrustAcceleration(t) - 9.80655)];
% end
% --------------------------------------------------------------------------

function [value,isterminal,direction] = events(t,y)
    % Locate the time when height passes through zero
    % and stop integration.
    value = y(1) > 0 || t < 3.45;     % detect height > 0, when height <= 0 gives 0
    isterminal = 1;   % stop the integration
    direction = 0;   % negative direction
end