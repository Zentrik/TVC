function [K_x, k_0, t, S] = affine_controller(state, current_u, horizon)
  persistent Ari Bri q w Inertia Mb
  
  if isempty(Ari)
    q = sym('q', [4 1], 'real');
    qi = sym('qi', [3 1], 'real');
    w = sym('w', [3 1], 'real');
    Mb = sym('Mb', [3 1]);
    I = sym('I', [3 3]);     

    qdot = .5 * quatmultiply(q.', [0;w].').';

    derivative_vector = [I \ (Mb - cross(w,I*w)); qdot(2:4); q(2:4)];

    Ari = jacobian(derivative_vector, [w; q(2:4); qi]);
    Bri = jacobian(derivative_vector, Mb);

    Inertia = diag([0.0829427911790996, 0.0829427911790996, 0.000246795169917015]);

    Ari = subs(Ari, I, Inertia);
    Bri = subs(Bri, I, Inertia);
  end

  quat = [(1 - norm(state(4:6))^2)^0.5; state(4:6)];
  ang = state(1:3);
  A = subs(Ari, q, quat);
  A = subs(A, w, ang);
  B = subs(Bri, q, quat);
  A = subs(Ari, Mb, current_u);

  A = double(A);
  B = double(B);

  qdot = .5 * quatmultiply(quat', [0, ang'])';
  c = [Inertia \ (current_u - cross(ang, Inertia*ang)); qdot(2:4); quat(2:4)];
  
  Q = diag([1,1,1, 2,2,2, 1,1,1]);
  R = diag([0.3 0.3 0.6]);
  [K_x, k_0, t, S] = Affine_LQR(A, B, Q, R, zeros(9, 1), current_u, state, current_u, horizon, c); % target current u, as we want to minise change in u.
end