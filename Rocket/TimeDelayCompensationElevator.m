R = 0.09022556390977443;
K_v = 47.29869062302208;
K_t = 0.01819548872180451;
m = 6.803886;
r = 0.02762679089;
G = 42.0 / 12.0 * 40.0 / 14.0;
A = [0 1; 1 -G^2 * K_t / (R*r^2*m*K_v)];
B = [0; G * K_t / (R*r*m)];
C = eye(2);
D = [0; 0];

L = 0.05;
dt = 0.005;
end_time = 10;
openloop_elevator = ss(A, B, C, D, dt);

discrete = expm([A, B; zeros(width(B), width(A) + width(B))] * dt);
Ad = discrete(1:length(A), 1:width(A));
Bd = discrete(1:length(B), width(A)+1:end);
K = dlqr(Ad, Bd, diag(1 ./ [0.02, 0.4].^2), 1/ 12^2);
Kt = K * (Ad - Bd*K)^(L/dt);
Kff = pinv(Bd);

% continuous feedforward and optimal tracking comparison (optimal tracking
% is better). L basically ignores Q and R so it has lower steady state
% error but has input which leads to higher overall cost. lim F as R goes
% to 0 is probably L.
% R = 0.01;
% Q = eye(2);
% [K, S] = lqr(A, B, Q, R);
% F = -inv(R) * B' * inv(A' - S*B*inv(R)*B') * Q;
% K + pinv(B)*A;
% L = [K 1] * pinv([A B; C D])*[zeros(2, 2); eye(2)] ;

% %if current state is ref
% ref = [1.524; 0];
% A*ref + B * -inv(R) * B' * (S + s_x) * ref
% A*ref + B * (-K * (ref - ref) + pinv(B)*-A*ref)

rem = [];
for t = 0:dt:end_time
  rem = [rem ref(t)];
end
plot(0:dt:end_time, rem(1, :))
hold on

x = [0; 0]; % initial state
for t = 0:dt:L - dt
%   K0 = K * (Ad - Bd*K)^(t/dt) * Ad^(L/dt - t/dt);
  u = K * (Ad - Bd*K)^(t/dt) * (ref(t) - x(:, end)) + Kff * (ref(t + dt) - Ad * ref(t));
  if abs(u) > 12
    u = 12 * sign(u);
    endy
  x = [Ad*x(:, 1) + Bd*u, x];
end
for t = L:dt:end_time - dt
  y = x(:, L / dt + 1);
  u = Kt * (ref(t) - y) + Kff * (ref(t + dt) - Ad * ref(t));
  if abs(u) > 12
    u = 12 * sign(u);
  end
  x = [Ad*x(:, 1) + Bd*u, x];
  %x = x(:, 1:L/sample_time + 1);
end
x = flip(x,2);
plot(0:dt:end_time, x(1, :))

function r = ref(t)
  if t < 5.1 && t > 0.05
    r = [1.524; 0];
  else
    r = [0;0];
  end
end 
