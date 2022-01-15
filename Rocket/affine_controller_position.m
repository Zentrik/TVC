function [K_x, k_0, t, S] = affine_controller_position(state, current_u, horizon, Thrust_magnitude)
  % for debugging affine_controller_position([rand(3, 1); zeros(3, 1); rand(4, 1)], rand(3, 1), inf, 20)
  
  persistent A_p B_p q w Inertia Mb Thrust Ve Xe
  
  if isempty(A_p)
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

    Aexyz = g + q2dcm(q(1), -q(2), -q(3), -q(4)) * ([Thrustxy; Thrustz]) / mass; % has to be inverse quat, g doesn't actually affect the devative_vector
    qdot = .5 * quatmultiply(q.', [0;w].').';

    devative_vector = [I \ (Mb - cross(w,I*w)); qdot(2:4) ; Aexyz(1:2); Ve]; % devative of x
    x = [w; q(2:4); Ve; Xe];
    u = Mb;

    A_p = jacobian(devative_vector, x);
    B_p = jacobian(devative_vector, u);

    Inertia = diag([0.0826975856, 0.0826975856, 2.4778e-04]);

    A_p = subs(A_p, I, Inertia);
    B_p = subs(B_p, I, Inertia);
  end

  quat = [(1 - norm(state(4:6))^2)^0.5; state(4:6)];
  ang = state(1:3);
  vel = state(7:8);

  A = subs(A_p, q, quat);
  A = subs(A, w, ang);
  A = subs(A, Thrust, Thrust_magnitude);
  A = subs(A, Mb, current_u);

  B = subs(B_p, q, quat);
  B = subs(B, Mb, current_u);
  B = subs(B, Thrust, Thrust_magnitude);
  
  A = double(A);
  B = double(B);

  qdot = .5 * quatmultiply(quat', [0, ang'])';

  Average_cg = 0.549556499983870; %trapz(masscgI(:, 1), masscgI(:, 5)) / masscgI(end, 1)
rocket_length = 0.95;
  r = [0; 0; -(rocket_length - Average_cg)]; % position of motor, point where thrust force oginates
  mass = 1.061675354603240; % trapz(masscgI(:, 1), masscgI(:, 2)) / masscgI(end, 1)
  g = [0; 0; -9.80655];
  
  Thrustxy = [0 1; -1 0] * current_u(1:2) / r(3);
  Thrustz = sqrt(Thrust_magnitude^2 - norm(Thrustxy)^2);

  Aexyz = g + quat2dcm(quatconj(quat')) * ([Thrustxy; Thrustz]) / mass;

  c = [Inertia \ (current_u - cross(ang, Inertia*ang)); qdot(2:4); Aexyz(1:2); vel];
  
  Q = diag([1,1,1, 2,2,2, 1,1, 1,1]);
  R = diag([0.3 0.3 0.47]);

  desired_state = zeros(10, 1);
  desired_control = current_u; % % target current u, as we want to minise change in u.
  [K_x, k_0, t, S] = Affine_LQR(A, B, Q, R, desired_state, desired_control, state, current_u, horizon, c);
end