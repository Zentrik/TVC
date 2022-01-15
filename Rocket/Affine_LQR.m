function [K_x, k_0, t, S_mat] = Affine_LQR(A, B, Q, R, options)
  % u = -Kx * x - k_0
  [A, B, Q, R, opt] = InputSantiser(A, B, Q, R, options);
  [nu, nx] = sizes(A, B);
  
  %A = 1; B = 1; c = 0; Q = 1; R = 1; Qf = Q; N = 0; t_final = 10;
  
  if opt.t_span_ode(1) == inf
    Q_xx = Q(1:nx, 1:nx);
    R_uu = R(1:nu, 1:nu);
    r_u = R(1:nu, nu + 1:end);
    q_x = Q(1:nx, nx + 1:end);
    
    %[S_xx, K_x] = icare(A, B, Q_xx, R_uu, N);
    [K_x, S_xx] = lqr(A, B, Q_xx, R_uu, opt.N);
    s_x = (-(N + S_xx * B) / R_uu * B' + A') \ ((N + S_xx * B) \ R_uu * r_u - S_xx * opt.c - q_x);
    % q_0 = Q(nx + 1:end, nx + 1:end);
    % r_0 = R(nu + 1:end, nu + 1:end);
    % q_0 + r_0 - (r_u + B' * s_x)' / R_uu * (r_u + B' * s_x) + 2 * s_x' * c
    S = [S_xx s_x; s_x' 0];
    [~, k_0] = affineKsoln(A, S, R, B, opt.N);
    S_mat = S;
    t = inf;
  else
    [t, S] = ode45(@(t, S) affineSdynamics(t, A, B, opt.c, Q, R, S, opt.N), opt.t_span_ode, opt.Qf);
    %S = reshape(S, nx + 1, nx + 1, length(t));
    %S = cellfun(@(v) reshape(v, nx + 1, nx + 1),mat2cell(S, ones(1,length(S)), (nx + 1)^2),'Uniform', 0);

    L = length(t);
    K_x = zeros(nu, nx, L);
    k_0 = zeros(nu, 1, L);
    S_mat = zeros(nx + 1, nx +1, L);
    for k = 1:length(t)
      S_mat(:, :, k) = reshape(S(k, :), nx + 1, nx + 1);
      [K_x(:, :, k), k_0(:, :, k)] = affineKsoln(t(k), A, S_mat(:, :, k), R, B, opt.N);
    end

    S_mat = flip(S_mat, 3);
    K_x = flip(K_x, 3);
    k_0 = flip(k_0, 3);
    t = flip(t);

%     plot(t, reshape(K_x(1, 4, :), 1, L))
  end
  
  for k = 1:length(t)
    k_0(:, :, k) = k_0(:, :, k) - opt.u_0(t(k)) - K_x(:, :, k) * opt.x_0(t(k));
  end
  
%   k_0 = k_0 - opt.u_0 - pagemtimes(K_x, opt.x_0);
end

function Sdot = affineSdynamics(t, A, B, c, Q, R, S, N)  
  [nu, nx] = sizes(A, B);

  S = reshape(S, nx + 1, nx + 1);
  
  S_xx = S(1:nx, 1:nx);
  s_x  = S(1:nx, nx + 1:end);
  s_0  = S(nx + 1:end, nx + 1:end);
  
  Q_xx = @(t) INDEX(Q(t), 1:nx, 1:nx);
  q_x  = @(t) INDEX(Q(t), 1:nx, nx + 1);
  q_0  = @(t) INDEX(Q(t), nx + 1, nx + 1);

  R_uu = @(t) INDEX(R(t), 1:nu, 1:nu);
  r_u  = @(t) INDEX(R(t), 1:nu, nu + 1);
  r_0  = @(t) INDEX(R(t), nu + 1, nu + 1);
  Ri   = @(t) inv(R_uu(t));
  
  if (min(eig(S_xx))<0) 
    warning('Drake:TVLQR:NegativeS','S is not positive definite'); 
  end
  
  Sdot_xx = -(Q_xx(t) - (N(t) + S_xx*B(t)) * Ri(t) * (N(t) + S_xx*B(t))' + S_xx*A(t) + A(t)'*S_xx);
  rs = @(t) (r_u(t) + B(t)'*s_x);
  Sdot_x = -(q_x(t) - (N(t) + S_xx*B(t)) * Ri(t) * rs(t) + A(t)'*s_x + S_xx*c(t));
  Sdot_0 = -(q_0(t) + r_0(t) - rs(t)'*Ri(t)*rs(t) + 2 * s_x'*c(t)); % I think you can just set this to 0 no problem
    
  Sdot = [Sdot_xx, Sdot_x; Sdot_x', Sdot_0];
  Sdot = reshape(Sdot, (nx + 1)^2, 1);
end

function [K_x, k_0] = affineKsoln(t, A, S, R, B, N)
  [nu, nx] = sizes(A, B);

  R_uu = @(t) INDEX(R(t), 1:nu, 1:nu);
  r_u  = @(t) INDEX(R(t), 1:nu, nu + 1);
  
  S_xx = S(1:nx, 1:nx);
  s_x  = S(1:nx, nx + 1:end);
  
  K_x = R_uu(t) \ (B(t)' * S_xx + N(t)');
  k_0 = R_uu(t) \ (B(t)' * s_x + r_u(t));
end

function [nu, nx] = sizes(A, B)
  nx = length(A(0)); nu = size(B(0), 2);
end

function [A, B, Q, R, options] = InputSantiser(A, B, Q, R, options)
  if ~isa(A, 'function_handle') % i.e. A not a function, so just a matrix
    A = @(t) A;
  end
  if ~isa(B, 'function_handle')
    B = @(t) B;
  end
  if ~isa(Q, 'function_handle')
    Q = @(t) Q;
  end
  if ~isa(R, 'function_handle')
    R = @(t) R;
  end
  
  [nu, nx] = sizes(A, B);
  
  if ~exist('options.x_0', 'var')
    options.x_0 = zeros(nx, 1);
  end
  if ~exist('options.u_0', 'var')
    options.u_0 = zeros(nu, 1);
  end
  if ~exist('options.x_d', 'var')
    options.x_d = options.x_0;
  end
  if ~exist('options.u_d', 'var')
    options.u_d = options.u_0;
  end
  if ~exist('options.c', 'var')
    options.c = zeros(nx, 1);
  end
  if ~exist('options.N', 'var')
    options.N = zeros(nx, nu);
  end
  
  if ~isa(options.x_0, 'function_handle')
    options.x_0 = @(t) options.x_0;
  end
  if ~isa(options.u_0, 'function_handle')
    options.u_0 = @(t) options.u_0;
  end
  if ~isa(options.x_d, 'function_handle')
    options.x_d = @(t) options.x_d;
  end
  if ~isa(options.u_d, 'function_handle')
    options.u_d = @(t) options.u_d;
  end
  if ~isa(options.c, 'function_handle')
    options.c = @(t) options.c;
  end
  if ~isa(options.N, 'function_handle')
    options.N = @(t) options.N;
  end
  
  if ~exist('options.Qf', 'var')
    options.Qf = zeros(nx+1);
  end
  
  [Q, R, options] = check(A, B, Q, R, options);
  
  options.x_d = @(t) options.x_d(t) - options.x_0(t); % desired state in linearised relative cooridnate system
  options.u_d = @(t) options.u_d(t) - options.u_0(t);
  
  Q_xx = @(t) INDEX(Q(t), 1:nx, 1:nx);
  q_x  = @(t) INDEX(Q(t), 1:nx, nx + 1);
  q_0  = @(t) INDEX(Q(t), nx + 1, nx + 1);

  R_uu = @(t) INDEX(R(t), 1:nu, 1:nu);
  r_u  = @(t) INDEX(R(t), 1:nu, nu + 1);
  r_0  = @(t) INDEX(R(t), nu + 1, nu + 1);
  
  q_x_tmp = @(t) q_x(t) - Q_xx(t) * options.x_d(t) - options.N(t) * options.u_d(t);
  q_0_tmp = @(t) q_0(t) + options.x_d(t)' * Q_xx(t) * options.x_d(t) + 2 * options.x_d(t)' * options.N(t) * options.u_d(t);
  
  Q = @(t) [Q_xx(t) q_x_tmp(t); q_x_tmp(t)' q_0_tmp(t)];
  
  r_u_tmp = @(t) r_u(t) - R_uu(t) * options.u_d(t) - options.N(t)' * options.x_d(t);
  r_0_tmp = @(t) r_0(t) + options.u_d(t)' * R_uu(t) * options.u_d(t);
  
  R = @(t) [R_uu(t) r_u_tmp(t); r_u_tmp(t)' r_0_tmp(t)];
end

function [Q, R, options] = check(A, B, Q, R, options)
  [nu, nx] = sizes(A, B);

  error(abcdchk(A(0), B(0)));

  if isa(Q(0),'double')
    if isequal(size(Q(0)), [nx nx])
      Q = @(t) [Q(t) zeros(nx, 1); zeros(1, nx) 0];
    end
    sizecheck(Q(0),[nx + 1, nx + 1], 'Q');
  else
    error('Q must be a double');
  end
  
  if isa(R(0),'double')
    if isequal(size(R(0)), [nu nu])
      R = @(t) [R(t) zeros(nu, 1); zeros(1, nu) 0];
    end
    sizecheck(R(0),[nu + 1, nu + 1], 'R');
  else
    error('R must be a double');
  end

  if isa(options.Qf,'double')
    sizecheck(options.Qf,[nx + 1,nx + 1], 'Qf');
  else
    error('Qf must be a double');
  end

  if isa(options.N(0),'double')
    sizecheck(options.N(0),[nx,nu], 'N');
  else
    error('N must be a double');
  end
  
  if isa(options.c(0),'double')
    if isequal(size(options.c(0)), [1 nx])
      options.c = @(t) options.c(t)';
    end
    sizecheck(options.c(0), [nx,1], 'options.c');
  else
    error('c must be a double');
  end
  
  if isa(options.x_d(0),'double')
    if isequal(size(options.x_d(0)), [1 nx])
      options.x_d = @(t) options.x_d(t)';
    end
    sizecheck(options.x_d(0), [nx,1], 'options.x_d');
  else
    error('x_d must be a double');
  end
  
  if isa(options.u_d(0),'double')
    if isequal(size(options.u_d(0)), [1 nu])
      options.u_d = @(t) options.u_d(t)';
    end
    sizecheck(options.u_d(0), [nu,1], 'options.u_d');
  else
    error('u_d must be a double');
  end
  
  if isa(options.x_0(0),'double')
    if isequal(size(options.x_0(0)), [1 nx])
      options.x_0 = @(t) options.x_0(t)';
    end
    sizecheck(options.x_0(0), [nx,1], 'options.x_0');
  else
    error('x_0 must be a double');
  end
  
  if isa(options.u_0(0),'double')
    if isequal(size(options.u_0(0)), [1 nu])
      options.u_0 = @(t)options.u_0(t)';
    end
    sizecheck(options.u_0(0), [nu,1], 'options.u_0');
  else
    error('u_0 must be a double');
  end
end

function sizecheck(Variable, Size, Name)
  if ~isequal(size(Variable), Size)
    error(join([Name, 'should be', '[', string(Size(1)), ',', string(Size(2)), ']', 'Size']));
  end
end

function out = INDEX(Matrix, R, C)
  out =  Matrix(R,C);
end