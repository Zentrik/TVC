function dxdt = RocketDerivative(x, u)
Mb = u(1:3);
Thrust = u(4);
I = diag(u(5:7));
COTCOM = u(8);

w = x(1:3);
q = [(1 - norm(x(4:6))^2)^0.5; x(4:6)];
Ve = x(7:9);
Xe = x(10:12);

dxdt = [I \ (Mb - cross(w,I*w)); .5 * quatmultiply(q.', [0;w].').'; quat2dcm([q(1), -q(2), -q(3), -q(4)]) * ([0;0;Thrust] + [Mb(1:2);0] / COTCOM) + [0;0;-9.80655]; Ve];
