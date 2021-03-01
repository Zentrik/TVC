function [K_x, k_0, S] = affine_controller(state)
  persistent Ari Bri q w Inertia
  
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
    Ari = subs(Ari, Mb, [0;0;0]);

    Bri = subs(Bri, I, Inertia);
  end

  quat = [(1 - norm(state(4:6))^2)^0.5; state(4:6)];
  ang = state(1:3);
  A = subs(Ari, q, quat);
  A = subs(A, w, ang);
  B = subs(Bri, q, quat);

  A = double(A);
  B = double(B);

  qdot = .5 * quatmultiply(quat', [0, ang'])';
  c = [Inertia \ ([0;0;0] - cross(ang,Inertia*ang)); qdot(2:4); quat(2:4)];
  
  Q = diag([1,1,1, 2,2,2, 1,1,1]);
  R = diag([0.3 0.3 200]);
  [K_x, k_0, ~, S] = Affine_LQR(A, B, Q, R, zeros(9, 1), zeros(3, 1), state, zeros(3, 1), inf, c);
end