function x = vdpStateFcnTranslationOnly(x, w, u) 
dt = 0.001; % [s] Sample time
x = x + vdpStateFcnContinuous(x, w, u) * dt;
if x(6) < 0
    x(3) = 0;
    x(6) = 0;
end
end

function dxdt = vdpStateFcnContinuous(x, wk, u)
Mb = u(1:3);
ThrustAcceleration = u(4);
COTCOM = u(5);
q = u(6:9);

Ve = x(1:3);
Xe = x(4:6);

qin = q';
DCM = [qin(:,1).^2 + qin(:,2).^2 - qin(:,3).^2 - qin(:,4).^2 , 2.*(qin(:,2).*qin(:,3) + qin(:,1).*qin(:,4)), 2.*(qin(:,2).*qin(:,4) - qin(:,1).*qin(:,3)) ; ...
    2.*(qin(:,2).*qin(:,3) - qin(:,1).*qin(:,4)), qin(:,1).^2 - qin(:,2).^2 + qin(:,3).^2 - qin(:,4).^2, 2.*(qin(:,3).*qin(:,4) + qin(:,1).*qin(:,2)) ; ...
    2.*(qin(:,2).*qin(:,4) + qin(:,1).*qin(:,3)), 2.*(qin(:,3).*qin(:,4) - qin(:,1).*qin(:,2)), qin(:,1).^2 - qin(:,2).^2 - qin(:,3).^2 + qin(:,4).^2];

dxdt = [DCM' * ([0;0;ThrustAcceleration] + [Mb(1:2);0] / COTCOM) + [0;0;-9.80655] + wk; Ve];
    
end