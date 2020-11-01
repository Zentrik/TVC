global Thrust
ThrustTable = xlsread('F15_Thrust','Sheet1');
Thrust = griddedInterpolant([-1; ThrustTable(:,1); 100], [0; ThrustTable(:,2); 0]);

numv = 5;
v = linspace(0, -100,numv);
nums = 2;
s = linspace(0.1, 100, nums);

data = zeros(numv, nums);
rocket(s(1), v(1))

% for j = 1:numv
%     for i = 1:nums
%         data(v(j), s(i)) = rocket(s(i), v(j));
%     end
% end

% surf(v, s, data)
% xlabel('Initial Velocity (m/s)');
% ylabel('Initial Height (m)');
% colormap(flipud(winter))

function impactVelocity = rocket(s0, v0)
tstart = 0;
tfinal = 30;
% plot(-2:0.001:6, Thrust(-2:0.001:6))
y0 = [s0; v0];
options = odeset('Events',@events,'Refine',4, 'RelTol', 1e-5, 'AbsTol', 1e-8);

tout = tstart;
yout = y0.';
teout = [];
yeout = [];
ieout = [];

t = 0;
nt = 1;
ye_t = 0;

% Solve until the first terminal event.
[t,y,te,ye,ie] = ode45(@f,[tstart tfinal],y0,options)

% Accumulate output.  This could be passed out as output arguments.
nt = length(t);
tout = [tout; t(2:nt)];
yout = [yout; y(2:nt,:)];
teout = [teout; te];          % Events at tstart are never reported.
yeout = [yeout; ye];
ieout = [ieout; ie];

% Set the new initial conditions, with .9 attenuation.
y0 = [0;0];
tstart = t(nt);
plot(t,y)
end

function dydt = f(t,y)
global Thrust
dydt = [y(2); Thrust(t)/1.0565 - 9.80655];
end
% --------------------------------------------------------------------------

function [value,isterminal,direction] = events(t,y)
% Locate the time when height passes through zero in a decreasing direction
% and stop integration.
value = y(1) > 0;     % detect height > 0, when height <= 0 gives 0
isterminal = 1;   % stop the integration
direction = 0;   % negative direction
end