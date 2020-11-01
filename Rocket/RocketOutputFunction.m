function y = RocketOutputFunction(x,u)
Mb = u(1:3);
Thrust = u(4);
I = diag(u(5:7));

w = x(1:3);
q = [1 - norm(x(4:7))^2; x(4:6)];

A = [(I(1,3)*I(2,2)^2*w(2) - I(1,2)*I(2,3)^2*w(3) + I(1,3)*I(3,2)^2*w(2) - I(1,2)*I(3,3)^2*w(3) + I(1,1)*I(1,2)*I(2,3)*w(2) - I(1,1)*I(1,3)*I(2,2)*w(2) - 2*I(1,2)*I(2,1)*I(2,3)*w(1) + 2*I(1,3)*I(2,1)*I(2,2)*w(1) + I(1,1)*I(1,2)*I(3,3)*w(3) - I(1,1)*I(1,3)*I(3,2)*w(3) - I(1,2)*I(2,2)*I(2,3)*w(2) + I(1,3)*I(2,2)*I(2,3)*w(3) - 2*I(1,2)*I(3,1)*I(3,3)*w(1) + 2*I(1,3)*I(3,1)*I(3,2)*w(1) - I(1,2)*I(3,2)*I(3,3)*w(2) + I(2,1)*I(2,2)*I(3,3)*w(3) - I(2,1)*I(2,3)*I(3,2)*w(3) + I(1,3)*I(3,2)*I(3,3)*w(3) - I(2,2)*I(3,1)*I(3,3)*w(2) + I(2,3)*I(3,1)*I(3,2)*w(2))/(I(1,1)*I(2,2)*I(3,3) - I(1,1)*I(2,3)*I(3,2) - I(1,2)*I(2,1)*I(3,3) + I(1,2)*I(2,3)*I(3,1) + I(1,3)*I(2,1)*I(3,2) - I(1,3)*I(2,2)*I(3,1)), (w(3)*I(1,2)^2*I(3,3) + 2*I(2,3)*w(2)*I(1,2)^2 - 2*w(2)*I(1,2)*I(1,3)*I(2,2) - w(3)*I(1,2)*I(1,3)*I(3,2) + I(2,3)*w(3)*I(1,2)*I(1,3) - I(2,3)*w(1)*I(1,2)*I(2,2) - w(1)*I(1,2)*I(3,2)*I(3,3) + I(1,1)*I(2,3)*w(1)*I(1,2) - w(3)*I(1,3)^2*I(2,2) + w(1)*I(1,3)*I(2,2)^2 - I(1,1)*w(1)*I(1,3)*I(2,2) + w(1)*I(1,3)*I(3,2)^2 + w(3)*I(2,2)^2*I(3,3) - 2*w(2)*I(2,2)*I(3,2)*I(3,3) - I(2,3)*w(3)*I(2,2)*I(3,2) - w(3)*I(2,2)*I(3,3)^2 - I(3,1)*w(1)*I(2,2)*I(3,3) + 2*I(2,3)*w(2)*I(3,2)^2 + I(2,3)*w(3)*I(3,2)*I(3,3) + I(2,3)*I(3,1)*w(1)*I(3,2))/(I(1,1)*I(2,2)*I(3,3) - I(1,1)*I(2,3)*I(3,2) - I(1,2)*I(2,1)*I(3,3) + I(1,2)*I(2,3)*I(3,1) + I(1,3)*I(2,1)*I(3,2) - I(1,3)*I(2,2)*I(3,1)), -(- w(2)*I(1,2)^2*I(3,3) - w(2)*I(1,2)*I(1,3)*I(2,3) - 2*w(3)*I(1,2)*I(1,3)*I(3,3) + I(3,2)*w(2)*I(1,2)*I(1,3) + w(1)*I(1,2)*I(2,3)^2 + w(1)*I(1,2)*I(3,3)^2 - I(1,1)*w(1)*I(1,2)*I(3,3) + w(2)*I(1,3)^2*I(2,2) + 2*I(3,2)*w(3)*I(1,3)^2 - w(1)*I(1,3)*I(2,2)*I(2,3) - I(3,2)*w(1)*I(1,3)*I(3,3) + I(1,1)*I(3,2)*w(1)*I(1,3) - w(2)*I(2,2)^2*I(3,3) - 2*w(3)*I(2,2)*I(2,3)*I(3,3) + I(3,2)*w(2)*I(2,2)*I(2,3) + w(2)*I(2,2)*I(3,3)^2 - I(2,1)*w(1)*I(2,2)*I(3,3) + 2*I(3,2)*w(3)*I(2,3)^2 - I(3,2)*w(2)*I(2,3)*I(3,3) + I(2,1)*I(3,2)*w(1)*I(2,3))/(I(1,1)*I(2,2)*I(3,3) - I(1,1)*I(2,3)*I(3,2) - I(1,2)*I(2,1)*I(3,3) + I(1,2)*I(2,3)*I(3,1) + I(1,3)*I(2,1)*I(3,2) - I(1,3)*I(2,2)*I(3,1)), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
 -(w(2)*I(1,1)^2*I(2,3) + w(3)*I(1,1)^2*I(3,3) - 2*w(1)*I(1,1)*I(2,1)*I(2,3) - I(1,3)*w(2)*I(1,1)*I(2,1) - w(3)*I(1,1)*I(2,3)^2 - I(2,2)*w(2)*I(1,1)*I(2,3) - 2*w(1)*I(1,1)*I(3,1)*I(3,3) - I(1,3)*w(3)*I(1,1)*I(3,1) - w(3)*I(1,1)*I(3,3)^2 - I(3,2)*w(2)*I(1,1)*I(3,3) + w(3)*I(2,1)^2*I(3,3) + 2*I(1,3)*w(1)*I(2,1)^2 - w(3)*I(2,1)*I(2,3)*I(3,1) + I(1,3)*w(3)*I(2,1)*I(2,3) - w(2)*I(2,1)*I(3,1)*I(3,3) + I(1,3)*I(2,2)*w(2)*I(2,1) + w(2)*I(2,3)*I(3,1)^2 + 2*I(1,3)*w(1)*I(3,1)^2 + I(1,3)*w(3)*I(3,1)*I(3,3) + I(1,3)*I(3,2)*w(2)*I(3,1))/(I(1,1)*I(2,2)*I(3,3) - I(1,1)*I(2,3)*I(3,2) - I(1,2)*I(2,1)*I(3,3) + I(1,2)*I(2,3)*I(3,1) + I(1,3)*I(2,1)*I(3,2) - I(1,3)*I(2,2)*I(3,1)), -(I(1,1)^2*I(2,3)*w(1) - I(1,3)^2*I(2,1)*w(3) + I(2,3)*I(3,1)^2*w(1) - I(2,1)*I(3,3)^2*w(3) - I(1,1)*I(1,3)*I(2,1)*w(1) + 2*I(1,1)*I(1,2)*I(2,3)*w(2) - 2*I(1,2)*I(1,3)*I(2,1)*w(2) + I(1,1)*I(1,3)*I(2,3)*w(3) - I(1,1)*I(2,2)*I(2,3)*w(1) + I(1,3)*I(2,1)*I(2,2)*w(1) + I(1,1)*I(1,2)*I(3,3)*w(3) - I(1,2)*I(1,3)*I(3,1)*w(3) - I(1,1)*I(3,2)*I(3,3)*w(1) + I(1,3)*I(3,1)*I(3,2)*w(1) + I(2,1)*I(2,2)*I(3,3)*w(3) - I(2,2)*I(2,3)*I(3,1)*w(3) - I(2,1)*I(3,1)*I(3,3)*w(1) - 2*I(2,1)*I(3,2)*I(3,3)*w(2) + 2*I(2,3)*I(3,1)*I(3,2)*w(2) + I(2,3)*I(3,1)*I(3,3)*w(3))/(I(1,1)*I(2,2)*I(3,3) - I(1,1)*I(2,3)*I(3,2) - I(1,2)*I(2,1)*I(3,3) + I(1,2)*I(2,3)*I(3,1) + I(1,3)*I(2,1)*I(3,2) - I(1,3)*I(2,2)*I(3,1)), (- w(1)*I(1,1)^2*I(3,3) - w(2)*I(1,1)*I(1,3)*I(2,3) - 2*w(3)*I(1,1)*I(1,3)*I(3,3) + I(3,1)*w(1)*I(1,1)*I(1,3) + w(1)*I(1,1)*I(2,3)^2 + w(1)*I(1,1)*I(3,3)^2 - I(1,2)*w(2)*I(1,1)*I(3,3) + w(2)*I(1,3)^2*I(2,1) + 2*I(3,1)*w(3)*I(1,3)^2 - w(1)*I(1,3)*I(2,1)*I(2,3) - I(3,1)*w(1)*I(1,3)*I(3,3) + I(1,2)*I(3,1)*w(2)*I(1,3) - w(1)*I(2,1)^2*I(3,3) - 2*w(3)*I(2,1)*I(2,3)*I(3,3) + I(3,1)*w(1)*I(2,1)*I(2,3) + w(2)*I(2,1)*I(3,3)^2 - I(2,2)*w(2)*I(2,1)*I(3,3) + 2*I(3,1)*w(3)*I(2,3)^2 - I(3,1)*w(2)*I(2,3)*I(3,3) + I(2,2)*I(3,1)*w(2)*I(2,3))/(I(1,1)*I(2,2)*I(3,3) - I(1,1)*I(2,3)*I(3,2) - I(1,2)*I(2,1)*I(3,3) + I(1,2)*I(2,3)*I(3,1) + I(1,3)*I(2,1)*I(3,2) - I(1,3)*I(2,2)*I(3,1)), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
 (w(2)*I(1,1)^2*I(2,2) + w(3)*I(1,1)^2*I(3,2) - 2*w(1)*I(1,1)*I(2,1)*I(2,2) - I(1,2)*w(2)*I(1,1)*I(2,1) - w(2)*I(1,1)*I(2,2)^2 - I(2,3)*w(3)*I(1,1)*I(2,2) - 2*w(1)*I(1,1)*I(3,1)*I(3,2) - I(1,2)*w(3)*I(1,1)*I(3,1) - w(2)*I(1,1)*I(3,2)^2 - I(3,3)*w(3)*I(1,1)*I(3,2) + w(3)*I(2,1)^2*I(3,2) + 2*I(1,2)*w(1)*I(2,1)^2 - w(3)*I(2,1)*I(2,2)*I(3,1) + I(1,2)*w(2)*I(2,1)*I(2,2) - w(2)*I(2,1)*I(3,1)*I(3,2) + I(1,2)*I(2,3)*w(3)*I(2,1) + w(2)*I(2,2)*I(3,1)^2 + 2*I(1,2)*w(1)*I(3,1)^2 + I(1,2)*w(2)*I(3,1)*I(3,2) + I(1,2)*I(3,3)*w(3)*I(3,1))/(I(1,1)*I(2,2)*I(3,3) - I(1,1)*I(2,3)*I(3,2) - I(1,2)*I(2,1)*I(3,3) + I(1,2)*I(2,3)*I(3,1) + I(1,3)*I(2,1)*I(3,2) - I(1,3)*I(2,2)*I(3,1)), -(- w(1)*I(1,1)^2*I(2,2) - 2*w(2)*I(1,1)*I(1,2)*I(2,2) - w(3)*I(1,1)*I(1,2)*I(3,2) + I(2,1)*w(1)*I(1,1)*I(1,2) + w(1)*I(1,1)*I(2,2)^2 - I(1,3)*w(3)*I(1,1)*I(2,2) + w(1)*I(1,1)*I(3,2)^2 + w(3)*I(1,2)^2*I(3,1) + 2*I(2,1)*w(2)*I(1,2)^2 - I(2,1)*w(1)*I(1,2)*I(2,2) - w(1)*I(1,2)*I(3,1)*I(3,2) + I(1,3)*I(2,1)*w(3)*I(1,2) + w(3)*I(2,2)^2*I(3,1) - w(1)*I(2,2)*I(3,1)^2 - 2*w(2)*I(2,2)*I(3,1)*I(3,2) - I(3,3)*w(3)*I(2,2)*I(3,1) - I(2,1)*w(3)*I(2,2)*I(3,2) + I(2,1)*w(1)*I(3,1)*I(3,2) + 2*I(2,1)*w(2)*I(3,2)^2 + I(2,1)*I(3,3)*w(3)*I(3,2))/(I(1,1)*I(2,2)*I(3,3) - I(1,1)*I(2,3)*I(3,2) - I(1,2)*I(2,1)*I(3,3) + I(1,2)*I(2,3)*I(3,1) + I(1,3)*I(2,1)*I(3,2) - I(1,3)*I(2,2)*I(3,1)), (I(1,1)^2*I(3,2)*w(1) - I(1,2)^2*I(3,1)*w(2) + I(2,1)^2*I(3,2)*w(1) - I(2,2)^2*I(3,1)*w(2) + I(1,1)*I(1,3)*I(2,2)*w(2) - I(1,2)*I(1,3)*I(2,1)*w(2) - I(1,1)*I(1,2)*I(3,1)*w(1) + I(1,1)*I(1,2)*I(3,2)*w(2) - I(1,1)*I(2,2)*I(2,3)*w(1) + I(1,2)*I(2,1)*I(2,3)*w(1) + 2*I(1,1)*I(1,3)*I(3,2)*w(3) - 2*I(1,2)*I(1,3)*I(3,1)*w(3) - I(2,1)*I(2,2)*I(3,1)*w(1) - I(1,1)*I(3,2)*I(3,3)*w(1) + I(1,2)*I(3,1)*I(3,3)*w(1) + I(2,1)*I(2,2)*I(3,2)*w(2) + 2*I(2,1)*I(2,3)*I(3,2)*w(3) - 2*I(2,2)*I(2,3)*I(3,1)*w(3) - I(2,1)*I(3,2)*I(3,3)*w(2) + I(2,2)*I(3,1)*I(3,3)*w(2))/(I(1,1)*I(2,2)*I(3,3) - I(1,1)*I(2,3)*I(3,2) - I(1,2)*I(2,1)*I(3,3) + I(1,2)*I(2,3)*I(3,1) + I(1,3)*I(2,1)*I(3,2) - I(1,3)*I(2,2)*I(3,1)), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
 q(1)/2, -q(4)/2, q(3)/2, 0, w(3)/2, -w(2)/2, 0, 0, 0, 0, 0, 0, 0;
 q(4)/2, q(1)/2, -q(2)/2, -w(3)/2, 0, w(1)/2, 0, 0, 0, 0, 0, 0, 0;
 -q(3)/2, q(2)/2, q(1)/2, w(2)/2, -w(1)/2, 0, 0, 0, 0, 0, 0, 0, 0;
 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0;
 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0;
 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0;
 0, 0, 0, (18014398509481984*Mb(1)*q(2))/3614073432593295 + (18014398509481984*Mb(2)*q(3))/3614073432593295 + 2*Thrust*q(4), (18014398509481984*Mb(2)*q(2))/3614073432593295 - (18014398509481984*Mb(1)*q(3))/3614073432593295 + 2*Thrust*q(1), 2*Thrust*q(2) - (18014398509481984*Mb(1)*q(4))/3614073432593295 - (18014398509481984*Mb(2)*q(1))/3614073432593295, 0, 0, 0, 0, 0, 0, 0;
 0, 0, 0, (18014398509481984*Mb(1)*q(3))/3614073432593295 - (18014398509481984*Mb(2)*q(2))/3614073432593295 - 2*Thrust*q(1), (18014398509481984*Mb(1)*q(2))/3614073432593295 + (18014398509481984*Mb(2)*q(3))/3614073432593295 + 2*Thrust*q(4), (18014398509481984*Mb(1)*q(1))/3614073432593295 - (18014398509481984*Mb(2)*q(4))/3614073432593295 + 2*Thrust*q(3), 0, 0, 0, 0, 0, 0, 0;
 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0;
 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0];

B = [(I(2,2)*I(3,3) - I(2,3)*I(3,2))/(I(1,1)*I(2,2)*I(3,3) - I(1,1)*I(2,3)*I(3,2) - I(1,2)*I(2,1)*I(3,3) + I(1,2)*I(2,3)*I(3,1) + I(1,3)*I(2,1)*I(3,2) - I(1,3)*I(2,2)*I(3,1)), -(I(1,2)*I(3,3) - I(1,3)*I(3,2))/(I(1,1)*I(2,2)*I(3,3) - I(1,1)*I(2,3)*I(3,2) - I(1,2)*I(2,1)*I(3,3) + I(1,2)*I(2,3)*I(3,1) + I(1,3)*I(2,1)*I(3,2) - I(1,3)*I(2,2)*I(3,1)), (I(1,2)*I(2,3) - I(1,3)*I(2,2))/(I(1,1)*I(2,2)*I(3,3) - I(1,1)*I(2,3)*I(3,2) - I(1,2)*I(2,1)*I(3,3) + I(1,2)*I(2,3)*I(3,1) + I(1,3)*I(2,1)*I(3,2) - I(1,3)*I(2,2)*I(3,1));
 -(I(2,1)*I(3,3) - I(2,3)*I(3,1))/(I(1,1)*I(2,2)*I(3,3) - I(1,1)*I(2,3)*I(3,2) - I(1,2)*I(2,1)*I(3,3) + I(1,2)*I(2,3)*I(3,1) + I(1,3)*I(2,1)*I(3,2) - I(1,3)*I(2,2)*I(3,1)), (I(1,1)*I(3,3) - I(1,3)*I(3,1))/(I(1,1)*I(2,2)*I(3,3) - I(1,1)*I(2,3)*I(3,2) - I(1,2)*I(2,1)*I(3,3) + I(1,2)*I(2,3)*I(3,1) + I(1,3)*I(2,1)*I(3,2) - I(1,3)*I(2,2)*I(3,1)), -(I(1,1)*I(2,3) - I(1,3)*I(2,1))/(I(1,1)*I(2,2)*I(3,3) - I(1,1)*I(2,3)*I(3,2) - I(1,2)*I(2,1)*I(3,3) + I(1,2)*I(2,3)*I(3,1) + I(1,3)*I(2,1)*I(3,2) - I(1,3)*I(2,2)*I(3,1));
 (I(2,1)*I(3,2) - I(2,2)*I(3,1))/(I(1,1)*I(2,2)*I(3,3) - I(1,1)*I(2,3)*I(3,2) - I(1,2)*I(2,1)*I(3,3) + I(1,2)*I(2,3)*I(3,1) + I(1,3)*I(2,1)*I(3,2) - I(1,3)*I(2,2)*I(3,1)), -(I(1,1)*I(3,2) - I(1,2)*I(3,1))/(I(1,1)*I(2,2)*I(3,3) - I(1,1)*I(2,3)*I(3,2) - I(1,2)*I(2,1)*I(3,3) + I(1,2)*I(2,3)*I(3,1) + I(1,3)*I(2,1)*I(3,2) - I(1,3)*I(2,2)*I(3,1)), (I(1,1)*I(2,2) - I(1,2)*I(2,1))/(I(1,1)*I(2,2)*I(3,3) - I(1,1)*I(2,3)*I(3,2) - I(1,2)*I(2,1)*I(3,3) + I(1,2)*I(2,3)*I(3,1) + I(1,3)*I(2,1)*I(3,2) - I(1,3)*I(2,2)*I(3,1));
 0, 0, 0;
 0, 0, 0;
 0, 0, 0;
 0, 0, 0;
 0, 0, 0;
 0, 0, 0;
 (9007199254740992*q(1)^2)/3614073432593295 + (9007199254740992*q(2)^2)/3614073432593295 - (9007199254740992*q(3)^2)/3614073432593295 - (9007199254740992*q(4)^2)/3614073432593295, (18014398509481984*q(2)*q(3))/3614073432593295 - (18014398509481984*q(1)*q(4))/3614073432593295, 0;
 (18014398509481984*q(1)*q(4))/3614073432593295 + (18014398509481984*q(2)*q(3))/3614073432593295, (9007199254740992*q(1)^2)/3614073432593295 + (9007199254740992*q(3)^2)/3614073432593295 - (9007199254740992*q(2)^2)/3614073432593295 - (9007199254740992*q(4)^2)/3614073432593295, 0;
 0, 0, 0;
 0, 0, 0];

C = eye(size(A,1));
D = zeros(size(B));

sys = absorbDelay(ss(A, B, C, D, 0.02, 'InputDelay',[3 3 1]));

y = sys.C * x + sys.D * Mb;