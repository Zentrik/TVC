function [K_x, k_0, t, S_mat] = Affine_LQR(A, B, Q, R, varargin)
  % u = -Kx * x - k_0
  [Q, R, c, t_final, x_d, u_d, x_0, u_0, N, Qf, nx, nu] = InputSantiser(A, B, Q, R, varargin);
  
  %A = 1; B = 1; c = 0; Q = 1; R = 1; Qf = Q; N = 0; t_final = 10;
  
  if t_final == inf
    Q_xx = Q(1:nx, 1:nx);
    R_uu = R(1:nu, 1:nu);
    r_u = R(1:nu, nu + 1:end);
    q_x = Q(1:nx, nx + 1:end);
    
    %[S_xx, K_x] = icare(A, B, Q_xx, R_uu, N);
    [K_x, S_xx] = lqr(A, B, Q_xx, R_uu, N);
    s_x = (-(N + S_xx * B) / R_uu * B' + A') \ ((N + S_xx * B) * inv(R_uu) * r_u - S_xx * c - q_x);
    % q_0 = Q(nx + 1:end, nx + 1:end);
    % r_0 = R(nu + 1:end, nu + 1:end);
    % q_0 + r_0 - (r_u + B' * s_x)' / R_uu * (r_u + B' * s_x) + 2 * s_x' * c
    S = [S_xx s_x; s_x' 0];
    [~, k_0] = affineKsoln(A, S, R, B, N);
    S_mat = S;
    t = inf;
  else
    [t, S] = ode45(@(t, S) affineSdynamics(A, B, c, Q, R, S, N), [t_final 0], Qf);
    %S = reshape(S, nx + 1, nx + 1, length(t));
    %S = cellfun(@(v) reshape(v, nx + 1, nx + 1),mat2cell(S, ones(1,length(S)), (nx + 1)^2),'Uniform', 0);

    L = length(t);
    K_x = zeros(nu, nx, L);
    k_0 = zeros(nu, 1, L);
    S_mat = zeros(nx + 1, nx +1, L);
    for i = 1:length(t)
      S_mat(:, :, i) = reshape(S(i, :), nx + 1, nx + 1);
      [K_x(:, :, i), k_0(:, :, i)] = affineKsoln(A, S_mat(:, :, i), R, B, N);
    end

    S_mat = flip(S_mat, 3);
    K_x = flip(K_x, 3);
    k_0 = flip(k_0, 3);
    t = flip(t);

    plot(t, reshape(K_x(1, 4, :), 1, L))
  end
  
  k_0 = k_0 - u_0 - pagemtimes(K_x, x_0);
end

function Sdot = affineSdynamics(A, B, c, Q, R, S, N)  
  [nu, nx] = sizes(A, B);

  S = reshape(S, nx + 1, nx + 1);
  
  S_xx = S(1:nx, 1:nx);
  s_x = S(1:nx, nx + 1:end);
  s_0 = S(nx + 1:end, nx + 1:end);
  
  Q_xx = Q(1:nx, 1:nx);
  q_x = Q(1:nx, nx + 1:end);
  q_0 = Q(nx + 1:end, nx + 1:end);

  R_uu = R(1:nu, 1:nu);
  r_u = R(1:nu, nu + 1:end);
  r_0 = R(nu + 1:end, nu + 1:end);
  Ri = inv(R_uu);
  
  if (min(eig(S_xx))<0) 
    warning('Drake:TVLQR:NegativeS','S is not positive definite'); 
  end
  
  Sdot_xx = -(Q_xx - (N+S_xx*B)*Ri*(N+S_xx*B)' + S_xx*A + A'*S_xx);
  rs = (r_u+B'*s_x);
  Sdot_x = -(q_x - (N+S_xx*B)*Ri*rs + A'*s_x + S_xx*c);
  Sdot_0 = -(q_0 + r_0 - rs'*Ri*rs + 2 * s_x'*c); % I think you can just set this to 0 no problem
    
  Sdot = [Sdot_xx, Sdot_x; Sdot_x', Sdot_0];
  Sdot = reshape(Sdot, (nx + 1)^2, 1);
end

function [K_x, k_0] = affineKsoln(A, S, R, B, N)
  [nu, nx] = sizes(A, B);

  R_uu = R(1:nu, 1:nu);
  r_u = R(1:nu, nu + 1:end);
  
  S_xx = S(1:nx, 1:nx);
  s_x = S(1:nx, nx + 1:end);
  
  K_x = R_uu \ (B' * S_xx + N');
  k_0 = R_uu \ (B' * s_x + r_u);
end

function [nu, nx] = sizes(A, B)
  nx = length(A); nu = size(B, 2);
end

function [Q, R, c, t_final, x_d, u_d, x_0, u_0, N, Qf, nx, nu] = InputSantiser(A, B, Q, R, extra_inputs)
  input_size = size(extra_inputs);
  if input_size(1) == 1
    if input_size(2) >= 1
      x_d = extra_inputs{1};
    end
    if input_size(2) >= 2
      u_d = extra_inputs{2};
    end
    if input_size(2) >= 3
      x_0 = extra_inputs{3};
    end
    if input_size(2) >= 4
      u_0 = extra_inputs{4};
    end
    if input_size(2) >= 5
      t_final = extra_inputs{5};
    end
    if input_size(2) >= 6
      c = extra_inputs{6};
    end
     if input_size(2) >= 7
      N = extra_inputs{7};
    end
    if input_size(2) >= 8
      Qf = extra_inputs{8};
    end
  end
  
  [nu, nx] = sizes(A, B);
  if ~exist('c', 'var')
    c = zeros(nx, 1);
  end
  if ~exist('t_final', 'var')
    t_final = inf;
  end
  if ~exist('x_d', 'var')
    x_d = zeros(nx, 1);
  end
  if ~exist('u_d', 'var')
    u_d = zeros(nu, 1);
  end
  if ~exist('x_0', 'var')
    x_0 = zeros(nx, 1);
  end
  if ~exist('u_0', 'var')
    u_0 = zeros(nu, 1);
  end
  if ~exist('N', 'var')
    N = zeros(nx, nu);
  end
  if ~exist('Qf', 'var')
    Qf = zeros(nx+1);
  end
  [Q, R, c, x_d, u_d, x_0, u_0] = check(A, B, c, Q, R, Qf, N, t_final, x_d, u_d, x_0, u_0);
  
  x_d = x_d - x_0; % desired state in linearised relative cooridnate system
  u_d = u_d - u_0;
  
  Q_xx = Q(1:nx, 1:nx);
  q_x = Q(1:nx, nx + 1:end);
  q_0 = Q(nx + 1:end, nx + 1:end);

  R_uu = R(1:nu, 1:nu);
  r_u = R(1:nu, nu + 1:end);
  r_0 = R(nu + 1:end, nu + 1:end);
  
  Q(1:nx, nx + 1:end) = q_x - Q_xx * x_d - N * u_d;
  Q(nx + 1:end, 1:nx) = Q(1:nx, nx + 1:end)';
  Q(nx + 1:end, nx + 1:end) = q_0 + x_d' * Q_xx * x_d + 2 * x_d' * N * u_d;
  
  R(1:nu, nu + 1:end) = r_u - R_uu * u_d - N' * x_d;
  R(nu + 1:end, 1:nu) = R(1:nu, nu + 1:end)';
  R(nu + 1:end, nu + 1:end) = r_0 + u_d' * R_uu * u_d;
end

function [Q, R, c, x_d, u_d, x_0, u_0] = check(A, B, c, Q, R, Qf, N, t_final, x_d, u_d, x_0, u_0)
  [nu, nx] = sizes(A, B);

  error(abcdchk(A, B));

  if isa(Q,'double')
    if size(Q) == [nx nx]
      Q = [Q zeros(nx, 1); zeros(1, nx) 0];
    end
    sizecheck(Q,[nx + 1,nx + 1], 'Q');
  else
    error('Q must be a double');
  end

  if (isa(R,'double'))
    if size(R) == [nu nu]
      R = [R zeros(nu, 1); zeros(1, nu) 0];
    end
    sizecheck(R,[nu + 1,nu + 1], 'R');
  else
    error('R must be a double');
  end

  if isa(Qf,'double')
    sizecheck(Qf,[nx + 1,nx + 1], 'Qf');
  else
    error('Qf must be a double');
  end

  if isa(N,'double')
    sizecheck(N,[nx,nu], 'N');
  else
    error('N must be a double');
  end
  
  if isa(c,'double')
    if size(c) == [1 nx]
      c = c';
    end
    sizecheck(c,[nx,1], 'c');
  else
    error('c must be a double');
  end
  
  if isa(t_final,'double')
    sizecheck(t_final, [1,1], 't_final');
  else
    error('t_final must be a double');
  end
  
  if isa(x_d,'double')
    if size(x_d) == [1 nx]
      x_d = x_d';
    end
    sizecheck(x_d, [nx,1], 'x_d');
  else
    error('x_d must be a double');
  end
  
  if isa(u_d,'double')
    if size(u_d) == [1 nu]
      u_d = u_d';
    end
    sizecheck(u_d, [nu,1], 'u_d');
  else
    error('u_d must be a double');
  end
  
  if isa(x_0,'double')
    if size(x_0) == [1 nx]
      x_0 = x_0';
    end
    sizecheck(x_0, [nx,1], 'x_0');
  else
    error('x_0 must be a double');
  end
  
  if isa(u_0,'double')
    if size(u_0) == [1 nu]
      u_0 = u_0';
    end
    sizecheck(u_0, [nu,1], 'u_0');
  else
    error('u_0 must be a double');
  end
end

function sizecheck(Variable, Size, Name)
  if ~isequal(size(Variable), Size)
    error(join([Name, 'should be', '[', string(Size(1)), ',', string(Size(2)), ']', 'Size']));
  end
end