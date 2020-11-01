function x = vdpStateFcnAttitudeOnly(xprev, w, u) 
dt = 0.001; % [s] Sample time
x = xprev + vdpStateFcnContinuous(xprev, w, u) * dt;
x(4:7) = x(4:7) / norm(x(4:7));
end

function dxdt = vdpStateFcnContinuous(x, wk, u)
Mb = u(1:3);
I = diag(u(4:6));

w = x(1:3);
q = x(4:7);

dxdt = [I \ (Mb + wk(1:3) - cross(w,I*w)); .5 * quatmultiply(q', [0;w]')'];
end