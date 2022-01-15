function [K_x, k_0, t, S] = Tracker(horizon, ThrustTable, desired_state, desired_moments)
  % for debugging affine_controller_position([rand(3, 1); zeros(3, 1); rand(4, 1)], rand(3, 1), inf, 20)
  
  persistent A_tmp B_tmp
  
  if isempty(A_tmp)
    q = sym('q', [4 1], 'real');
    w = sym('w', [3 1], 'real');
    I = sym('I', [3 3]);     
    Ve = sym('Ve', [3 1]);
    Xe = sym('Xe', [3 1]);

    Mb = sym('Moment', [3 1], 'real');
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

    devative_vector = [I \ (Mb - cross(w,I*w)); qdot(2:4) ; Aexyz; Ve]; % devative of x
    x = [w; q(2:4); Ve; Xe];
    u = Mb;

    A_p = jacobian(devative_vector, x);
    B_p = jacobian(devative_vector, u);

    Inertia = diag([0.0826975856, 0.0826975856, 2.4778e-04]);

    A_p = subs(A_p, I, Inertia);
    B_p = subs(B_p, I, Inertia);
    
    A_tmp = matlabFunction(A_p, 'Vars', {[Xe; Ve; q; w], Mb, Thrust});
    B_tmp = matlabFunction(B_p, 'Vars', {[Xe; Ve; q; w], Mb, Thrust});
  end
  
  Q = diag([1,1,0.5, 2,2,1, 1,1,1, 1,1,1]);
  R = diag([1 1 10]);
  
%   u_d = @(time) interp1(desired_moments(:, 1), desired_moments(:, 2:end), time);
%   x_d = @(time) interp1(desired_state(:, 1), desired_state(:, 2:end), time);
%   Thrust_magnitude = @(time) max(interp1(ThrustTable(:, 1), ThrustTable(:, 2), time), eps); % (Thrust - |Moment|)^1/2 in B has problems if Thrust = 0

  initial_time = desired_moments(1, 1);
  step_size = desired_moments(2, 1) - initial_time;
  
  u_d = @(time) interpolate(desired_moments(:, 2:end), (time - initial_time) / step_size + 1);
  x_d = @(time) interpolate(desired_state(:, 2:end), (time - initial_time) / step_size + 1);
  
  initial_time_thrust = ThrustTable(1, 1);
  step_size_thrust = ThrustTable(2, 1) - initial_time_thrust;
  
  Thrust_magnitude = @(time) max(interpolate(ThrustTable(:, 2), (time - initial_time_thrust) / step_size_thrust + 1), 1e-2);

  A = @(t) A_tmp(x_d(t)', u_d(t)', Thrust_magnitude(t));
  B = @(t) B_tmp(x_d(t)', u_d(t)', Thrust_magnitude(t));
      
%   opt.x_0 = x_d;
%   opt.u_0 = u_d;
  opt.t_span_ode = [3.45 3.45 - horizon];
  opt.Qf = diag([10,10,5, 20,20,0, 100,100,100, 0,0, 100]);
  
%   [K_x, k_0, t, S] = Affine_LQR(A, B, Q, R, opt);
  [K_x, k_0, t, S] = Trajectory_LQR_fast(A, B, Q, R, opt); % k_0 should be 0 always, i.e. all(all(k_0 == 0)) == 1
%   for k = 1:length(t)
%     x = x_d(t(k));
%     k_0(:, :, k) = k_0(:, :, k) - u_d(t(k))' - K_x(:, :, k) * x([1:6, 8:13])';
%   end
end