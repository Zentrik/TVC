%% Attitude only
q = sym('q', [4 1], 'real');
qi = sym('gi', [3 1]);
w = sym('w', [3 1], 'real');
I = sym('I', [3 3]);     

Mb = sym('Moment', [3 1]);

qdot = .5 * quatmultiply(q.', [0;w].').';

derivative_vector = [I \ (Mb - cross(w,I*w)); qdot(2:4); q(2:4)]; % derivative of x
x = [w; q(2:4); qi];
u = Mb;

A = jacobian(derivative_vector, x);
B = jacobian(derivative_vector, u);
C = eye(size(A,1));
D = zeros(size(B));

Inertia = diag([0.0826975856, 0.0826975856, 2.4778e-04]);
omega = [0;0;0];

A = subs(A, I, Inertia);
A = subs(A, w, omega);
A = subs(A, q, [1; 0; 0; 0]);

B = subs(B, I, Inertia);
B = subs(B, q, [1; 0; 0; 0]);

A = double(A);
B = double(B);

Q = diag([1,1,1, 2,2,2, 1,1,1]);
R = diag([0.3 0.3 0.6]);
K = lqr(A, B, Q, R);
%K(abs(K)<1e-5)=0;

Kd = lqrd(A, B, Q, R, 0.02); % discrete lqr controller for continous system

roll_control_v= ss(A, B, C, D);
roll_control_discrete = c2d(roll_control_v, 0.02);
roll_control_delay = ss(A, B, C, D, 'InputDelay',[0.07 0.07 0.01]);
roll_control_delay_approx = minreal(pade(roll_control_delay, 3)); 
%[a,b,c,d,e] = ctrbf(roll_control_delay_approx.A, roll_control_delay_approx.B, roll_control_delay_approx.C); [sum(e), length(roll_control_delay_approx.A)]

Qy = roll_control_delay_approx.C' * Q * roll_control_delay_approx.C;
% lqr(roll_control_delay_approx, Qy, Rv), same as lqry
icare(roll_control_delay_approx.A, roll_control_delay_approx.B, Qy, R)
ry = lqry(roll_control_delay_approx, Q, R); % continous lqr controller based on output
Kyd = lqrd(roll_control_delay_approx.A, roll_control_delay_approx.B, Qy, R, 0.02 );
%observer_poles = real(eig(A - B * Krd)) * 5 + imag(eig(A - B * Krd))*1i;
%observer_poles =
%observer_poles(mod(0:length(roll_control_delay_approx.A)-1, numel(observer_poles)) + 1);

% Observer regulator tuned through pole placement
observer_poles = [-20:-1:-19 - length(roll_control_delay_approx.A)]';
L = place(roll_control_delay_approx.A', roll_control_delay_approx.C', observer_poles)';
controller = reg(roll_control_delay_approx, ry, L);   

% LQG controller
kalman_roll_control_delay_approx = minreal(pade(ss(A,[B B],C,0, 'InputDelay',[0.07 0.07 0.01 0 0 0]), 5));
[kest, L, P] = kalman(kalman_roll_control_delay_approx, diag([2 2 2/5]), 0.1);
kalman_controller = lqgreg(kest, ry);

Mb_bar = [0.1;0.1;0.01];
Vd = diag(roll_control_delay_approx.B * Mb_bar);
Vn = 0.01;
[L, P, E] = lqe(roll_control_delay_approx.A, roll_control_delay_approx.B, roll_control_delay_approx.C, diag(Mb_bar), Vn * eye(size(roll_control_delay_approx.C, 1)));
lqe_estimator = estim(roll_control_delay_approx, L, 1:13, 1:3);

lqr(roll_control_delay_approx.A', roll_control_delay_approx.C', Vd, Vn)';

%% BEST ESTIMATOR SO FAR
L_simulink = lqe(roll_control_delay_approx.A, eye(length(roll_control_delay_approx.A)), roll_control_delay_approx.C, 0.001 * eye(28), 0.0001 * eye(13));
simulink_estimator = estim(roll_control_delay_approx, L_simulink, 1:13, 1:3);

%{
for i = 1:10
    try
        kalman(kalman_roll_control_delay_approx, diag([i i i/5]), 0.1);
        i
    catch
    end
end
%}
%estim(roll_control_delay_approx, L, [1:13], [1:3])

%% No quaternion integral, velocity control
q = sym('q', [4 1], 'real');
w = sym('w', [3 1], 'real');
I = sym('I', [3 3]);     
Ve = sym('Ve', [2 1]);

Mb = sym('Moment', [3 1]);
Thrust = sym('Thrust', 'real');

syms q2dcm(qin);
qin = sym('qin', [1 4], 'real');

q2dcm(qin) = [qin(:,1).^2 + qin(:,2).^2 - qin(:,3).^2 - qin(:,4).^2 , 2.*(qin(:,2).*qin(:,3) + qin(:,1).*qin(:,4)), 2.*(qin(:,2).*qin(:,4) - qin(:,1).*qin(:,3)) ; 2.*(qin(:,2).*qin(:,3) - qin(:,1).*qin(:,4)), qin(:,1).^2 - qin(:,2).^2 + qin(:,3).^2 - qin(:,4).^2, 2.*(qin(:,3).*qin(:,4) + qin(:,1).*qin(:,2)) ; 2.*(qin(:,2).*qin(:,4) + qin(:,1).*qin(:,3)), 2.*(qin(:,3).*qin(:,4) - qin(:,1).*qin(:,2)), qin(:,1).^2 - qin(:,2).^2 - qin(:,3).^2 + qin(:,4).^2];

Average_cg = 0.549556499983870; %trapz(masscgI(:, 1), masscgI(:, 5)) / masscgI(end, 1)
rocket_length = 0.95;
r = [0; 0; -(rocket_length - Average_cg)]; % position of motor, point where thrust force oginates
mass = 1.061675354603240; % trapz(masscgI(:, 1), masscgI(:, 2)) / masscgI(end, 1)
g = [0; 0; -9.80655];

% r \times Thrust = Mb, if r = [0; 0; z] then Thrust = (Mb(2) / z, - Mb(1) / z, 0), then we add in z component of Thrust
% (Mb(2) / z, - Mb(1) / z) = [0 1; -1 0] * Mb / z

Thrustxy = [0 1; -1 0] * Mb(1:2) / r(3);
Thrustz = sqrt(Thrust^2 - norm(Thrustxy)^2);

Aexyz = g + q2dcm(q(1), -q(2), -q(3), -q(4)) * ([Thrustxy; Thrustz]) / mass; % has to be inverse quat, g doesn't actually affect the derivative_vector
qdot = .5 * quatmultiply(q.', [0;w].').';

derivative_vector = [I \ (Mb - cross(w,I*w)); qdot(2:4) ; Aexyz(1:2)]; % derivative of x
x = [w; q(2:4); Ve];
u = Mb;

A = jacobian(derivative_vector, x);
B = jacobian(derivative_vector, u);
C = eye(size(A,1));
D = zeros(size(B));

Inertia = diag([0.0826975856, 0.0826975856, 2.4778e-04]);
omega = [0;0;0];

A = subs(A, I, Inertia);
A = subs(A, w, omega);
A = subs(A, q, [1; 0; 0; 0]);
A = subs(A, Thrust, 10.6);
A = subs(A, Mb, [0; 0; 0]);

B = subs(B, I, Inertia);
B = subs(B, q, [1; 0; 0; 0]);

A = double(A);
B = double(B);

Q = diag([1,1,1, 2,2,2, 1,1]);
R = diag([0.3 0.3 0.47]);
K = lqr(A, B, Q, R);
%Krv(abs(Krv)<1e-5)=0;

Kd = lqrd(A, B, Q, R, 0.02); % discrete lqr controller for continous system

roll_control_v = ss(A, B, C, D);
roll_control_discrete = c2d(roll_control_v, 0.02);
roll_control_delay = ss(A, B, C, D, 'InputDelay',[0.07 0.07 0.01]);
roll_control_delay_approx = minreal(pade(roll_control_delay, 3)); 
%[a,b,c,d,e] = ctrbf(roll_control_delay_approx.A, roll_control_delay_approx.B, roll_control_delay_approx.C); [sum(e), length(roll_control_delay_approx.A)]

Qy = roll_control_delay_approx.C' * Q * roll_control_delay_approx.C;
% lqr(roll_control_delay_approx, Qy, Rv), same as lqry
%     icare(roll_control_delay_approx.A, roll_control_delay_approx.B, Qy, Rv)
ry = lqry(roll_control_delay_approx, Q, R); % continous lqr controller based on output
Kyd = lqrd(roll_control_delay_approx.A, roll_control_delay_approx.B, Qy, R, 0.02 );
%observer_poles = real(eig(A - B * Krd)) * 5 + imag(eig(A - B * Krd))*1i;
%observer_poles =
%observer_poles(mod(0:length(roll_control_delay_approx.A)-1, numel(observer_poles)) + 1);

% Observer regulator tuned through pole placement
observer_poles = [-20:-1:-19 - length(roll_control_delay_approx.A)]';
L = place(roll_control_delay_approx.A', roll_control_delay_approx.C', observer_poles)';
controller = reg(roll_control_delay_approx, ry, L);   

% LQG controllers
[kest, L, P] = kalman(roll_control_delay_approx, diag([2 2 2/5]), 0.1 * eye(size(roll_control_delay_approx.C, 1)));
kalman_controller = lqgreg(kest, ry);

Mb_bar = [0.1;0.1;0.01];
Vd = diag(roll_control_delay_approx.B * Mb_bar);
Vn = 0.01 * eye(size(roll_control_delay_approx.C, 1));
[L, P, E] = lqe(roll_control_delay_approx.A, roll_control_delay_approx.B, roll_control_delay_approx.C, diag(Mb_bar), Vn);
lqe_estimator = estim(roll_control_delay_approx, L, 1:size(roll_control_delay_approx.C, 1), 1:3);

% lqr(roll_control_delay_approx.A', roll_control_delay_approx.C', Vd, Vn)';  

%% No quaternion integral, velocity control with integral
q = sym('q', [4 1], 'real');
w = sym('w', [3 1], 'real');
I = sym('I', [3 3]);     
Ve = sym('Ve', [2 1]);
Xe = sym('Xe', [2 1]);

Mb = sym('Moment', [3 1]);
Thrust = sym('Thrust', 'real');

syms q2dcm(qin);
qin = sym('qin', [1 4], 'real');

q2dcm(qin) = [qin(:,1).^2 + qin(:,2).^2 - qin(:,3).^2 - qin(:,4).^2 , 2.*(qin(:,2).*qin(:,3) + qin(:,1).*qin(:,4)), 2.*(qin(:,2).*qin(:,4) - qin(:,1).*qin(:,3)) ; 2.*(qin(:,2).*qin(:,3) - qin(:,1).*qin(:,4)), qin(:,1).^2 - qin(:,2).^2 + qin(:,3).^2 - qin(:,4).^2, 2.*(qin(:,3).*qin(:,4) + qin(:,1).*qin(:,2)) ; 2.*(qin(:,2).*qin(:,4) + qin(:,1).*qin(:,3)), 2.*(qin(:,3).*qin(:,4) - qin(:,1).*qin(:,2)), qin(:,1).^2 - qin(:,2).^2 - qin(:,3).^2 + qin(:,4).^2];

Average_cg = 0.549556499983870; %trapz(masscgI(:, 1), masscgI(:, 5)) / masscgI(end, 1)
rocket_length = 0.95;
r = [0; 0; -(rocket_length - Average_cg)]; % position of motor, point where thrust force oginates
mass = 1.061675354603240; % trapz(masscgI(:, 1), masscgI(:, 2)) / masscgI(end, 1)
g = [0; 0; -9.80655];

% r \times Thrust = Mb, if r = [0; 0; z] then Thrust = (Mb(2) / z, - Mb(1) / z, 0), then we add in z component of Thrust
% (Mb(2) / z, - Mb(1) / z) = [0 1; -1 0] * Mb / z

Thrustxy = [0 1; -1 0] * Mb(1:2) / r(3);
Thrustz = sqrt(Thrust^2 - norm(Thrustxy)^2);

Aexyz = g + q2dcm(q(1), -q(2), -q(3), -q(4)) * ([Thrustxy; Thrustz]) / mass; % has to be inverse quat, g doesn't actually affect the derivative_vector
qdot = .5 * quatmultiply(q.', [0;w].').';

derivative_vector = [I \ (Mb - cross(w,I*w)); qdot(2:4) ; Aexyz(1:2); Ve]; % derivative of x
x = [w; q(2:4); Ve; Xe];
u = Mb;

% derivative_vector = [Ve; Aexyz(1:2); qdot(2:4); I \ (Mb - cross(w,I*w))]; % derivative of x
% x = [Xe; Ve; q(2:4); w];

A = jacobian(derivative_vector, x);
B = jacobian(derivative_vector, u);
C = eye(size(A,1));
D = zeros(size(B));

Inertia = diag([0.0826975856, 0.0826975856, 2.4778e-04]);
omega = [0;0;0];

A = subs(A, I, Inertia);
A = subs(A, w, omega);
A = subs(A, q, [1; 0; 0; 0]);
A = subs(A, Thrust, 10.6);
% A = subs(A, Thrust, 15.589230769230769);
A = subs(A, Mb, [0; 0; 0]);

B = subs(B, I, Inertia);
B = subs(B, q, [1; 0; 0; 0]);

A = double(A);
B = double(B);

Q = diag([1,1,1, 2,2,2, 1,1, 1,1]);
R = diag([0.3 0.3 0.47]);
K = lqr(A, B, Q, R);
%Krv(abs(Krv)<1e-5)=0;

Kd = lqrd(A, B, Q, R, 0.02); % discrete lqr controller for continous system

roll_control_v = ss(A, B, C, D);
roll_control_discrete = c2d(roll_control_v, 0.02);
roll_control_delay = ss(A, B, C, D, 'InputDelay',[0.07 0.07 0.01]);
roll_control_delay_approx = minreal(pade(roll_control_delay, 3)); 
%[a,b,c,d,e] = ctrbf(roll_control_delay_approx.A, roll_control_delay_approx.B, roll_control_delay_approx.C); [sum(e), length(roll_control_delay_approx.A)]

Qy = roll_control_delay_approx.C' * Q * roll_control_delay_approx.C;
% lqr(roll_control_delay_approx, Qy, Rv), same as lqry
%     icare(roll_control_delay_approx.A, roll_control_delay_approx.B, Qy, Rv)
ry = lqry(roll_control_delay_approx, Q, R); % continous lqr controller based on output
Kyd = lqrd(roll_control_delay_approx.A, roll_control_delay_approx.B, Qy, R, 0.02 );
%observer_poles = real(eig(A - B * Krd)) * 5 + imag(eig(A - B * Krd))*1i;
%observer_poles =
%observer_poles(mod(0:length(roll_control_delay_approx.A)-1, numel(observer_poles)) + 1);

% Observer regulator tuned through pole placement
observer_poles = [-20:-1:-19 - length(roll_control_delay_approx.A)]';
L = place(roll_control_delay_approx.A', roll_control_delay_approx.C', observer_poles)';
controller = reg(roll_control_delay_approx, ry, L);   

% LQG controllers
[kest, L, P] = kalman(roll_control_delay_approx, diag([2 2 2/5]), 0.1 * eye(size(roll_control_delay_approx.C, 1)));
kalman_controller = lqgreg(kest, ry);

Mb_bar = [0.1;0.1;0.01];
Vd = diag(roll_control_delay_approx.B * Mb_bar);
Vn = 0.01 * eye(size(roll_control_delay_approx.C, 1));
[L, P, E] = lqe(roll_control_delay_approx.A, roll_control_delay_approx.B, roll_control_delay_approx.C, diag(Mb_bar), Vn);
lqe_estimator = estim(roll_control_delay_approx, L, 1:size(roll_control_delay_approx.C, 1), 1:3);

% lqr(roll_control_delay_approx.A', roll_control_delay_approx.C', Vd, Vn)';  

%% No quaternion integral, velocity control with integral with HEIGHT
q = sym('q', [4 1], 'real');
w = sym('w', [3 1], 'real');
I = sym('I', [3 3]);     
Ve = sym('Ve', [3 1]);
Xe = sym('Xe', [3 1]);

Mb = sym('Moment', [3 1]);
Thrust = sym('Thrust', 'real');

syms q2dcm(qin);
qin = sym('qin', [1 4], 'real');

q2dcm(qin) = [qin(:,1).^2 + qin(:,2).^2 - qin(:,3).^2 - qin(:,4).^2 , 2.*(qin(:,2).*qin(:,3) + qin(:,1).*qin(:,4)), 2.*(qin(:,2).*qin(:,4) - qin(:,1).*qin(:,3)) ; 2.*(qin(:,2).*qin(:,3) - qin(:,1).*qin(:,4)), qin(:,1).^2 - qin(:,2).^2 + qin(:,3).^2 - qin(:,4).^2, 2.*(qin(:,3).*qin(:,4) + qin(:,1).*qin(:,2)) ; 2.*(qin(:,2).*qin(:,4) + qin(:,1).*qin(:,3)), 2.*(qin(:,3).*qin(:,4) - qin(:,1).*qin(:,2)), qin(:,1).^2 - qin(:,2).^2 - qin(:,3).^2 + qin(:,4).^2];

Average_cg = 0.549556499983870; %trapz(masscgI(:, 1), masscgI(:, 5)) / masscgI(end, 1)
rocket_length = 0.95;
r = [0; 0; -(rocket_length - Average_cg)]; % position of motor, point where thrust force oginates
mass = 1.061675354603240; % trapz(masscgI(:, 1), masscgI(:, 2)) / masscgI(end, 1)
g = [0; 0; -9.80655];

% r \times Thrust = Mb, if r = [0; 0; z] then Thrust = (Mb(2) / z, - Mb(1) / z, 0), then we add in z component of Thrust
% (Mb(2) / z, - Mb(1) / z) = [0 1; -1 0] * Mb / z

Thrustxy = [0 1; -1 0] * Mb(1:2) / r(3);
Thrustz = sqrt(Thrust^2 - norm(Thrustxy)^2);

Aexyz = g + q2dcm(q(1), -q(2), -q(3), -q(4)) * ([Thrustxy; Thrustz]) / mass; % has to be inverse quat, g doesn't actually affect the derivative_vector
qdot = .5 * quatmultiply(q.', [0;w].').';

derivative_vector = [I \ (Mb - cross(w,I*w)); qdot(2:4) ; Aexyz; Ve]; % derivative of x
x = [w; q(2:4); Ve; Xe];
u = Mb;

A = jacobian(derivative_vector, x);
B = jacobian(derivative_vector, u);
C = eye(size(A,1));
D = zeros(size(B));

Inertia = diag([0.0826975856, 0.0826975856, 2.4778e-04]);
omega = [0;0;0];

A = subs(A, I, Inertia);
A = subs(A, w, omega);
A = subs(A, q, [1; 0; 0; 0]);
A = subs(A, Thrust, 10.6);
A = subs(A, Mb, [0; 0; 0]);

B = subs(B, I, Inertia);
B = subs(B, q, [1; 0; 0; 0]);
B = subs(B, Mb, [0; 0; 0]);

A = double(A);
B = double(B);