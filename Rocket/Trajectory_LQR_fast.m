function [K_x, k_0, t, S_mat] = Trajectory_LQR_fast(A, B, Q, R, options)
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
    [t, S] = ode45(@(t, S) affineSdynamics(t, A, B, opt.c, Q, R, S, opt.N), opt.t_span_ode, reshape(opt.Qf, (nx + 1)^2, 1));
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
end

function Sdot = affineSdynamics(t, A, B, c, Q, R, S, N)  
  [nu, nx] = sizes(A, B);

  S = reshape(S, nx + 1, nx + 1);
  
  S_xx = S(1:nx, 1:nx);
  s_x  = S(1:nx, nx + 1:end);
  s_0  = S(nx + 1:end, nx + 1:end);
  
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
  
  Sdot_xx = -(Q_xx - (N + S_xx*B(t)) * Ri * (N + S_xx*B(t))' + S_xx*A(t) + A(t)'*S_xx);
  rs = @(t) (r_u + B(t)'*s_x);
  Sdot_x = -(q_x - (N + S_xx*B(t)) * Ri * rs(t) + A(t)'*s_x + S_xx*c);
  Sdot_0 = -(q_0 + r_0 - rs(t)'*Ri*rs(t) + 2 * s_x'*c); % I think you can just set this to 0 no problem
    
  Sdot = [Sdot_xx, Sdot_x; Sdot_x', Sdot_0];
  Sdot = reshape(Sdot, (nx + 1)^2, 1);
  
  if ~isreal(Sdot)
    warning('Sdot not real at ' + string(t))
  end
end

function [K_x, k_0] = affineKsoln(t, A, S, R, B, N)
  [nu, nx] = sizes(A, B);

  R_uu = R(1:nu, 1:nu);
  r_u  = R(1:nu, nu + 1:end);
  
  S_xx = S(1:nx, 1:nx);
  s_x  = S(1:nx, nx + 1:end);
  
  K_x = R_uu \ (B(t)' * S_xx' + N'); % R_uu \ (N + S_xx * B(t))'
  k_0 = R_uu \ (B(t)' * s_x + r_u);
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
  
  [nu, nx] = sizes(A, B);
  
  if ~isfield(options, 'x_0')
    options.x_0 = zeros(nx, 1);
  end
  if ~isfield(options, 'u_0')
    options.u_0 = zeros(nu, 1);
  end
  if ~isfield(options, 'x_d')
    options.x_d = options.x_0;
  end
  if ~isfield(options, 'u_d')
    options.u_d = options.u_0;
  end
  if ~isfield(options, 'c')
    options.c = zeros(nx, 1);
  end
  if ~isfield(options, 'N')
    options.N = zeros(nx, nu);
  end
  if ~isfield(options, 'Qf')
    options.Qf = zeros(nx+1);
  end
  
  [Q, R, options] = check(A, B, Q, R, options);
end

function [Q, R, options] = check(A, B, Q, R, options)
  [nu, nx] = sizes(A, B);

  error(abcdchk(A(0), B(0)));

  if isa(Q,'double')
    if isequal(size(Q), [nx nx])
      Q = [Q zeros(nx, 1); zeros(1, nx) 0];
    end
    sizecheck(Q,[nx + 1, nx + 1], 'Q');
  else
    error('Q must be a double');
  end
  
  if isa(R,'double')
    if isequal(size(R), [nu nu])
      R = [R zeros(nu, 1); zeros(1, nu) 0];
    end
    sizecheck(R,[nu + 1, nu + 1], 'R');
  else
    error('R must be a double');
  end

  if isa(options.Qf,'double')
    if isequal(size(options.Qf), [nx nx])
      options.Qf = [options.Qf zeros(nx, 1); zeros(1, nx) 0];
    end
    sizecheck(options.Qf,[nx + 1,nx + 1], 'Qf');
  else
    error('Qf must be a double');
  end

  if isa(options.N,'double')
    sizecheck(options.N,[nx,nu], 'N');
  else
    error('N must be a double');
  end
  
  if isa(options.c,'double')
    if isequal(size(options.c), [1 nx])
      options.c = options.c';
    end
    sizecheck(options.c, [nx,1], 'options.c');
  else
    error('c must be a double');
  end
  
  if isa(options.x_d,'double')
    if isequal(size(options.x_d), [1 nx])
      options.x_d = options.x_d';
    end
    sizecheck(options.x_d, [nx,1], 'options.x_d');
  else
    error('x_d must be a double');
  end
  
  if isa(options.u_d,'double')
    if isequal(size(options.u_d), [1 nu])
      options.u_d = options.u_d';
    end
    sizecheck(options.u_d, [nu,1], 'options.u_d');
  else
    error('u_d must be a double');
  end
  
  if isa(options.x_0,'double')
    if isequal(size(options.x_0), [1 nx])
      options.x_0 = options.x_0';
    end
    sizecheck(options.x_0, [nx,1], 'options.x_0');
  else
    error('x_0 must be a double');
  end
  
  if isa(options.u_0,'double')
    if isequal(size(options.u_0), [1 nu])
      options.u_0 = options.u_0';
    end
    sizecheck(options.u_0, [nu,1], 'options.u_0');
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