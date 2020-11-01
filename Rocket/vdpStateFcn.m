function x_next = vdpStateFcn(x, wk, u) 
dt = 0.001; % [s] Sample time

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

GlobalAcceleration = DCM' * Acceleration + [0;0;-9.80655];
if x(13) <= 0 && GlobalAcceleration(3) < 0
    GlobalAcceleration(3) = 0;
end

x_next(1:3) = w + dt * ( I \ (Mb + wk(1:3) - cross(w,I*w)));
x_next(4:7) = quatnormalize(quatmultiply(q, quatexp([0 dt/2 * (w + x_next(1:3))/2])));
x_next(8:10) = Ve + dt * GlobalAcceleration; 
x_next(11:13) = Xe + dt * Ve + dt^2/2 * GlobalAcceleration;

if x_next(13) < 0
    x_next(10) = 0;
    x_next(13) = 0;
end
end

% function dxdt = vdpStateFcnContinuous(x, wk, u, dt)
% Mb = u(1:3);
% Acceleration = u(4:6);
% I = diag(u(7:9));
% 
% w = x(1:3);
% q = x(4:7);
% Ve = x(8:10);
% Xe = x(11:13);
% gyrobias = x(14:16);
% accelbias = x(17:19);
% 
% qin = q';
% DCM = [qin(:,1).^2 + qin(:,2).^2 - qin(:,3).^2 - qin(:,4).^2 , 2.*(qin(:,2).*qin(:,3) + qin(:,1).*qin(:,4)), 2.*(qin(:,2).*qin(:,4) - qin(:,1).*qin(:,3)) ; ...
%     2.*(qin(:,2).*qin(:,3) - qin(:,1).*qin(:,4)), qin(:,1).^2 - qin(:,2).^2 + qin(:,3).^2 - qin(:,4).^2, 2.*(qin(:,3).*qin(:,4) + qin(:,1).*qin(:,2)) ; ...
%     2.*(qin(:,2).*qin(:,4) + qin(:,1).*qin(:,3)), 2.*(qin(:,3).*qin(:,4) - qin(:,1).*qin(:,2)), qin(:,1).^2 - qin(:,2).^2 - qin(:,3).^2 + qin(:,4).^2];
% 
% %dxdt = [I \ (Mb + wk(1:3) - cross(w,I*w)); .5 * quatmultiply(q', [0;w]')'; DCM' * ([0;0;Thrust] + [Mb(1:2);0] / COTCOM) / mass + [0;0;-9.80655]; Ve];
% 
% GlobalAcceleration = DCM' * Acceleration + [0;0;-9.80655];
% if x(13) <= 0 && GlobalAcceleration(3) < 0
%     GlobalAcceleration(3) = 0;
% end
% dxdt = [I \ (Mb + wk(1:3) - cross(w,I*w)); GlobalAcceleration; Ve + dt/2 * GlobalAcceleration; zeros(6,1)];
% % dxdt = [I \ (Mb + wk(1:3) - cross(w,I*w)); .5 * quatmultiply(q', [0;w]')'];
% end