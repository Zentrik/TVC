function [estimatedRollState, estimatedPitchState] = NonlinearStateEstimationForActuatorDelay(state, all_u, ThrustTable, time, delay, masscgI)
  %Inertia = diag([0.0829427911790996, 0.0829427911790996, 0.000246795169917015]);
  %{
  persistent Mb_buffer Mb_t 
  if isempty(Mb_buffer)
    Mb_t = [-0.5 -eps];
    Mb_buffer = zeros(3, 2); 
  end
  
  Mb_buffer = [Mb_buffer u];
  Mb_t = [Mb_t time];
  
  while time - Mb_t(2) > 0.5
    Mb_t = Mb_t(2:end);
    Mb_buffer = Mb_buffer(:, 2:end);
  end
  
  [Mb_t, idx]= unique(Mb_t);
  Mb_buffer = Mb_buffer(:, idx);
  
  Mb_t_corrected = Mb_t + 0.5 - time; %should be 0 for current actuator
  %}
  
  [~, roll_x] = ode45(@(t, x) derivative_vector(all_u, x, ThrustTable, time, t, masscgI, delay), [0 min(delay)], state);
  [~, pitch_x] = ode45(@(t, x) derivative_vector(all_u, x, ThrustTable, time, t, masscgI, delay), [min(delay) max(delay)], roll_x(end, :));
  estimatedRollState = roll_x(end, [1:3, 5:end]);
  estimatedPitchState = pitch_x(end, [1:3, 5:end]);
end

function x_dot = derivative_vector(all_u, x, ThrustTable, time, t, masscgI, delay)
  w = x(1:3);
  q = x(4:7);
  %qi = x(8:10);
  Ve = x(11:13);
  Xe = x(14:16);
  %Mb = [interp1(Mb_t, Mb_buffer(1, :)', t); interp1(Mb_t, Mb_buffer(2, :)', t); interp1(Mb_t, Mb_buffer(3, :)', t)];
  x = interpolate(all_u(:, 1), all_u(:, 2), time + t - delay(1));
  y = interpolate(all_u(:, 1), all_u(:, 3), time + t - delay(2));
  z = interpolate(all_u(:, 1), all_u(:, 4), time + t - delay(3));
  Mb = [x; y; z];
%   Mb = all_u(end, 2:4)';
%   Mb = interpolate(all_u(:, 1), all_u(:, 2:4), time + t)';
  
%   Ixx = interpolate(masscgI(:,1), masscgI(:,3), time + t);
%   Izz = interpolate(masscgI(:,1), masscgI(:,4), time + t);
  
  
  Thrust = interpolate(ThrustTable(:,1), ThrustTable(:, 2), time + t);
%   mass = interpolate(masscgI(:, 1), masscgI(:, 2), time + t);
%   cg = interpolate(masscgI(:, 1), masscgI(:, 5), time + t);
  
  temp = interpolate(masscgI(:, 1), masscgI(:, 2:5), time + t);
  mass = temp(1);
  Ixx = temp(2);
  Izz = temp(3);
  cg = temp(4);
  
  I = diag([Ixx, Ixx, Izz]);
  
  Mb(1:2) = normalise(Mb(1:2)) * min(Thrust * sin(deg2rad(5)), norm(Mb(1:2)) );
  if Xe(3) <= 0
    Mb = [0;0;0];
  end
  
  qdot = .5 * quatmultiply(q.', [0;w].').';
  
  Force = [0 1; -1 0] * Mb(1:2) / (0.95 - cg);
  Ae = quatrotate(quatinv(q'), [Force; (Thrust^2 - norm(Force)^2)^.5]')' / mass + [0;0;-9.7863774];
  if Xe(3) <= 0
    Ae(end-2:end-1) = [0;0];
    Ae(end) = max(0, Ae(end));
  end
  
  x_dot = [I \ (Mb - cross(w,I*w)); qdot; q(2:4); Ae; Ve];
end

function out = normalise(in)
  length = norm(in);
  if length ~= 0
    out = in / length;
  else
    out = in;
  end
end

function out = interpolate(x, y, xq)
  out = interp1(x, y, xq, 'linear');
  if sum(isnan(out)) >= 1
%     if xq > x(end)
%       out(isnan(out)) = y(end, isnan(out));
%     elseif xq < x(1)
%       out(isnan(out)) = y(1, isnan(out));
%     end
    if xq > x(end)
      out = y(end, :);
    elseif xq < x(1)
      out = y(1, :);
    end
  end
end