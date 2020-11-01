dt = 0.001;
I = diag([0.08 0.08 0.002]);
Mb = [0.01; 0.01; 0.0002];
w = [0; 0; 0];
P = 1e-4 * eye(6);
q = [1 0 0 0];

for t = 1:100
    alpha = I \ (Mb - cross(w,I*w));
    wmeas = x(1:3) + normrnd(0, [0.1; 0.1; 0.1], 3, 1);
    [q, P] = update(wmeas, alpha, q, P, dt) ;
end

function [q_t, P3] = update(y, alpha, x, P, dt) 
w = x(1:3);

% w_t = w + dt * alpha
w_t = w + dt * alpha;
n_t = n;
% q_t = q_t-1 * exp(dt/2 * (w_t-1 + w_t) /2)
q_t = quatnormalize(quatmultiply(q, quatexp([0 dt/2 * (w + w_t)/2'])));

F = eye(6);
G = eye(6);
Q = diag([0.001 0.001 0.00001 0.001 0.001 0.00001]);
%P = P_t - 1|t - 1
%P2 = P_t|t-1
P2 = F*P*F' + G*Q*G';

%epsilon = y_t - y_t|t-1, h(x_t|t-1) = w_t|t-1
epsilon = y - w_t;
H = [eye(3), zeros(3,3)];
R = diag([0.1 0.1 0.1]);
%P2 = P_t|t-1
S = H * P2 * H' + R;
K = P2 * H' * inv(S);
n = K*epsilon;

P3 = P2 - K * S * K';

q_t = quatmultiply(quatexp(n/2), q);
n = [0 0 0];
end