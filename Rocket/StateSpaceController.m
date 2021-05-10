function [nocontrol, polecontrol, control, LQRcontrol, LQIcontrol, noroll, norollpitch] = StateSpaceController(Inertia, omega, wbw)
    %% No Control
    % symbolic variables
    I = sym('I', [3 3]);
    %I = [I(1) 0 0;0 I(2,2) 0; 0 0 I(3,3)]; %diagonal I matrix
    InverseI = sym('InverseI', [3 3]);
    g = sym('q', [3 1], 'real');
    w = sym('w', [3 1], 'real');
    Mb = sym('Mb', [3 1]);
    
    % Define State Matrices
    A = jacobian([I \ (Mb - cross(w,I*w)); 0.5 * w], [w; g]);
    B = jacobian([I \ (Mb - cross(w,I*w)); 0.5 * w], Mb);
    C = eye(6);
    D = zeros([6 3]);
    
    A = subs(A, I, Inertia);
    A = subs(A, InverseI, inv(Inertia));
    
    B = subs(B, I, Inertia);
    B = subs(B, InverseI, inv(Inertia));
    
    A = subs(A, w, omega);
    B = subs(B, w, omega);
    
    %{
    A = subs(A, q, Quaternion);
    B = subs(B, q, Quaternion);
    %}
    
    A = double(A); 
    B = double(B); 
    C = double(C); 
    D = double(D);
    
    % Create State Space Object
    nocontrol = ss(A, B, C, D);
    
    % Open loop EigenValues
    eig(A);
    
    %% Desired Closed loop EigenValues, BAD Poles unstable system, high steady state erroris
    p = [0, 0, -5, -5, 0, -4].';
    K = place(A, B, p);
    
    Apcl = A - B*K;
    polecontrol = ss(Apcl, B, C, 0);
    
    %% Recommended Control from paper
    % Quaternion Controller
    Kp = Inertia * wbw^2;
    %Kp(:,3) = zeros(3,1); remove roll
    
    % angular velocity controller
    Kpd =  2*Inertia*wbw;
    %Kpd(:,3) = zeros(3,1);
    
    Kr = [Kpd Kp];
    Ky = [Kpd Kp];
    
    % Define Closed Loop State Matrices
    Acl = A - B*Ky*C;
    Bcl = B*Kr;
    Ccl = C;
    Dcl = zeros(size(C,1),size(Kr,2));
    
    % Define State Space Models
    control = ss(A - B*Ky*C, B*Kr, C, Dcl);
    
    %% LQR Control
    Q = diag([1, 1, 0, 2, 2, eps]);
    R = 1;
    [K, S, P] = lqr(nocontrol, Q, R);
    K(abs(K)<1e-5)=0;
    
    Alcl = A - B*K;
    LQRcontrol = ss(Alcl, B, C, 0);
    
    %% LQI control
    Ai = [A, zeros(size(A,1)); eye(size(A,1)),  zeros(size(A,1))];
    Bi = [B; zeros(size(A,1), size(B,2))];
    Ci = eye(size(Ai,1));
    Di = 0;
    
    LQIcontrol = ss(Ai, Bi, Ci, Di);
  
    Q = diag([.5,.5,eps,3,3,eps,.5,.5,eps,6,6,eps]);
    R = 1;
    
    Ki = lqr(LQIcontrol, Q, R);
    Ki(abs(Ki)<1e-5)=0;
    LQIcontrol = ss(Ai - Bi*Ki, Bi, Ci, Di);
        
    %% LQI, no w3 or q3 or w integral
    Inertia = 0.0826975856;
    A2 = [0 0 0 0;
          0 0 0 0;
          0.5 0 0 0;
          0 0.5 0 0];
    B2 = [1/Inertia(1) 0;
          0 1/Inertia(1);
          0 0;
          0 0];
    C2 = eye(4);
    D2 = 0;
    
    Ai2 = [A2, zeros(4,2);zeros(2), eye(2), zeros(2)];
    Bi2 = [B2; zeros(2)];
    Ci2 = eye(size(Ai2,1));
      
    Q = diag([1,1,5,5,3,3]);
    R = .01;
    Ki2 = lqr(Ai2, Bi2, Q, R);
    noroll = ss(Ai2 - Bi2*Ki2, Bi2, Ci2, D2);
    %%  No roll with Pitch moment
    syms Ixx F vinf m COT_TO_COM
    g = sym('q', [2 1]);
    gi = sym('gi', [2 1]);
    w = sym('w', [2 1]);
    Mb = sym('Mb', [2 1]);
    
    Reference_Diameter = 7.62e-2;
    Reference_Radius = Reference_Diameter / 2;
    Reference_Area = Reference_Radius^2 * pi;
    air_density = 1.2228;
    C_pitch = -1.3658;
    k = Reference_Area * Reference_Diameter * air_density * C_pitch;
    
    A2m = jacobian([(Mb + [k * (norm(Mb) + F)/m*vinf;0])/Ixx; 0.5 * w; g], [w; g;gi]); %one term to provide constant acceleration to free stream velocity
    A2m = subs(A2m, [Ixx F m], [Inertia(1), 0, 1.0567]);
    A2m = subs(A2m, Mb, [sqrt(0.5); sqrt(0.5)]); % so norm is 1, should be average after testing
    A2m = double(A2m);
    
    %approximate norm(Mb) to be sigma Mb
    B2m = jacobian([(Mb + [k * (Mb(1)/COT_TO_COM + Mb(2)/COT_TO_COM + F) / m * vinf;0])/Ixx; 0.5 * w; g], Mb);
    B2m = subs(B2m, [Ixx vinf m COT_TO_COM], [Inertia(1) 12 1.0567 .39] );
    B2m = double(B2m);
    C2m = eye(size(A2m,1));
    
    Q = diag([.5,0.5,1,1,.5,.5]);
    R = 1;
    Ki2m = lqr(A2m, B2m, Q, R);
    norollpitch = ss(A2m - B2m * Ki2m, B2m, C2m, 0); 
    
    %% no roll with dual integral
    syms Ixx
    g = sym('q', [2 1]);
    gi = sym('gi', [2 1]);
    w = sym('w', [2 1]);
    wi = sym('wi', [2 1]);
    Mb = sym('Mb', [2 1]);
    
    Aii2 = double(jacobian([Mb/ Ixx; 0.5 * w; w; g], [w; g; wi; gi]));
    Bii2 = jacobian([Mb/ Ixx; 0.5 * w; w; g], Mb);
    
    %% no roll, PIDish
    g = sym('q', [2 1]);
    gi = sym('gi', [2 1]);
    w = sym('w', [2 1]);
    Mb = sym('Mb', [2 1]);
    
    Ixx = 0.0826975856;
    Ai2 = double(jacobian([Mb/ Ixx; 0.5 * w; g], [w; g; gi]));
    Bi2 = double(jacobian([Mb/ Ixx; 0.5 * w; g], Mb));
    C2 = eye(size(Ai2,1));
    D2 = 0;
   
    Q = diag([1,1,5,5,3,3]);
    R = 0.1;
    Ki2 = lqr(Ai2, Bi2, Q, R);
    noroll = ss(Ai2 - Bi2*Ki2, Bi2, Ci2, D2);
    
    %% Precompensator tuner
    g = sym('q', [3 1]);
    w = sym('w', [3 1]);
    Mb = sym('Mb', [3 1]);
    I = sym('I', [3 3]);
    
    Inertia = diag([0.0826975856, 0.0826975856, 2.4778e-04]);
    omega = [0;0;0];
    
    Apre = jacobian([I \ (Mb - cross(w,I*w)); 0.5 * w], [w; g]);
    Bpre = jacobian([I \ (Mb - cross(w,I*w)); 0.5 * w], Mb);
    Cpre = eye(size(Apre,1));
    Dpre = zeros(size(Bpre));

    Apre = subs(Apre, I, Inertia);
    Bpre = subs(Bpre, I, Inertia);
    Apre = subs(Apre, w, omega);
    Bpre = subs(Bpre, w, omega);
    
    Apre = double(Apre);
    Bpre = double(Bpre);
    
    Q = diag([1, 1, .5, 2, 2, 1.5]);
    R = 1;
    Kpre = lqr(Apre, Bpre, Q, R);
    precomp = ss(Apre - Bpre*Kpre, Bpre, Cpre, Dpre);
    dcgain(precomp);
    step(precomp);
    N = 2^.5;
    precomp = ss(Apre - Bpre*Kpre, Bpre*N, Cpre, Dpre);
   %% horizontal velocity control
   Acc = sym('Acc', [2 1]);
   v = sym('v', [2 1]);
   x = sym('x', [2 1]);
   
   Av = double(jacobian([Acc; v], [v; x]));
   Bv = double(jacobian([Acc; v], Acc));
   C = eye(4);
   D = 0;
   
   Q = diag([1, 1, .5, .5]);
   R = 2;
   Kv = lqr(Av, Bv, Q, R);
   velocityss = ss(Av -Bv*Kv, Bv, C, D);
   
   %% Roll control integral
    g = sym('q', [3 1]);
    gi = sym('gi', [3 1]);
    w = sym('w', [3 1]);
    Mb = sym('Mb', [3 1]);
    I = sym('I', [3 3]);
    r = sym('r', [3 1]);
    
    %Inertia = diag([0.0826975856, 0.0826975856, 2.4778e-04]);
    Inertia = diag([0.0829427911790996, 0.0829427911790996, 0.000246795169917015]);
    omega = [0;0;0];
    
    Ari = jacobian([I \ (Mb - cross(w,I*w)); 0.5 * w; g - r], [w; g; gi]); %r is the reference quaternion vector
    Bri = jacobian([I \ (Mb - cross(w,I*w)); 0.5 * w; g - r], Mb);
    Bref = jacobian([I \ (Mb - cross(w,I*w)); 0.5 * w; g - r], r);
    Cri = eye(size(Ari,1));
    Dri = zeros(size(Bref));

    Ari = subs(Ari, I, Inertia);
    Bri = subs(Bri, I, Inertia);
    Ari = subs(Ari, w, omega);
    Bri = subs(Bri, w, omega);
    
    Ari = double(Ari);
    Bri = double(Bri);
    Bref = double(Bref);
    
    Q = diag([1, 1, 1, 5, 5, 5, 1, 1, 1]);
    R = 0.001;
    Kri = lqr(Ari, Bri, Q, R);
    roll_control_i= ss(Ari - Bri*Kri, Bref, Cri, Dri);
    
    %% roll control no integral
    g = sym('q', [3 1]);
    w = sym('w', [3 1]);
    Mb = sym('Mb', [3 1]);
    I = sym('I', [3 3]);
    
    Inertia = diag([0.0826975856, 0.0826975856, 2.4778e-04]);
    omega = [0;0;0];
    
    Ar = jacobian([I \ (Mb - cross(w,I*w)); 0.5 * w], [w; g]);
    Br = jacobian([I \ (Mb - cross(w,I*w)); 0.5 * w], Mb);
    Cr = eye(size(Ar,1));
    Dr = zeros(size(Br));

    Ar = subs(Ar, I, Inertia);
    Br = subs(Br, I, Inertia);
    Ar = subs(Ar, w, omega);
    Br = subs(Br, w, omega);
    
    Ar = double(Ar);
    Br = double(Br);
    
    Q = diag([1, 1, 1, 5, 5, 5]);
    R = .001;
    Kr = lqr(Ar, Br, Q, R);
    roll_control = ss(Ar - Br*Kr, Br, Cr, Dr);
    Kdc = dcgain(roll_control);
    N = 1/Kdc(4);
    roll_control =  ss(Ar - Br*Kr, Br * N, Cr, Dr);
    %% just roll
    Ajr = [0 0 0; 0.5 0 0; 0 1 0];
    Bjr = [1/2.4778e-04; 0; 0];
    Cjr = eye(size(Ajr,1));
    Djr = zeros(size(Bjr));
    
    Qjr = diag([1, 2, 1]);
    Rjr = .1;
    Kjr = lqr(Ajr, Bjr, Qjr, Rjr);
    
    just_roll = ss(Ajr, Bjr, Cjr, Djr);
    just_roll_control = ss(Ajr - Bjr * Kjr, Bjr, Cjr, Djr);
    
    %% Roll control integral with velocity
    q = sym('q', [4 1], 'real');
    qi = sym('gi', [3 1]);
    w = sym('w', [3 1], 'real');
    Mb = sym('Mb', [3 1]);
    I = sym('I', [3 3]);     
    Ve = sym('Ve', [2 1]);
    Xe = sym('Xe', [2 1]);
    Thrust = sym('Thrust', 'real');
    
    syms q2dcm(qin);
    qin = sym('qin', [1 4], 'real');
    
    q2dcm(qin) = [qin(:,1).^2 + qin(:,2).^2 - qin(:,3).^2 - qin(:,4).^2 , 2.*(qin(:,2).*qin(:,3) + qin(:,1).*qin(:,4)), 2.*(qin(:,2).*qin(:,4) - qin(:,1).*qin(:,3)) ; 2.*(qin(:,2).*qin(:,3) - qin(:,1).*qin(:,4)), qin(:,1).^2 - qin(:,2).^2 + qin(:,3).^2 - qin(:,4).^2, 2.*(qin(:,3).*qin(:,4) + qin(:,1).*qin(:,2)) ; 2.*(qin(:,2).*qin(:,4) + qin(:,1).*qin(:,3)), 2.*(qin(:,3).*qin(:,4) - qin(:,1).*qin(:,2)), qin(:,1).^2 - qin(:,2).^2 - qin(:,3).^2 + qin(:,4).^2];
    
    Aexyz =  q2dcm(q(1), -q(2), -q(3), -q(4)) * ([0;0;Thrust] + [Mb(1:2);0] / 0.401242753755115); % has to be inverse quat
    qdot = .5 * quatmultiply(q.', [0;w].').';
    
    derivative_vector = [I \ (Mb - cross(w,I*w)); qdot(2:4) ; q(2:4); Aexyz(1:2); Ve];
    
    Ari = jacobian(derivative_vector, [w; q(2:4); qi; Ve; Xe]);
    Bri = jacobian(derivative_vector, Mb);
    Cri = eye(size(Ari,1));
    Dri = zeros(size(Bri));
    
    Inertia = diag([0.0826975856, 0.0826975856, 2.4778e-04]);
    omega = [0;0;0];

    Ari = subs(Ari, I, Inertia);
    Ari = subs(Ari, w, omega);
    Ari = subs(Ari, q, [1; 0; 0; 0]);
    Ari = subs(Ari, Mb, [0;0;0]);
    Ari = subs(Ari, Thrust, 10.6);
    
    Bri = subs(Bri, I, Inertia);
    Bri = subs(Bri, q, [1; 0; 0; 0]);
    
    Ari = double(Ari);
    Bri = double(Bri);
    Bref = [Bri, [zeros(6,3);-eye(3);zeros(4,3)]];
    
    Qv = diag([1,1,1, 2,2,2, 1,1,1, 1,1, 1,1]);
    Rv = diag([0.3 0.3 0.6]);
    Krv = lqr(Ari, Bri, Qv, diag([0.3 0.3 0.47]) );
    %Krv(abs(Krv)<1e-5)=0;
    
    Krd = lqrd(Ari, Bri, Qv, Rv, 0.02); % discrete lqr controller for continous system
    
    roll_control_v= ss(Ari, Bri, Cri, Dri);
    roll_control_discrete = c2d(roll_control_v, 0.02);
    roll_control_delay = ss(Ari, Bri, Cri, Dri, 'InputDelay',[0.07 0.07 0.01]);
    roll_control_delay_approx = minreal(pade(roll_control_delay, 3)); 
    %[a,b,c,d,e] = ctrbf(roll_control_delay_approx.A, roll_control_delay_approx.B, roll_control_delay_approx.C); [sum(e), length(roll_control_delay_approx.A)]
    
    Qy = roll_control_delay_approx.C' * Qv * roll_control_delay_approx.C;
    % lqr(roll_control_delay_approx, Qy, Rv), same as lqry
    icare(roll_control_delay_approx.A, roll_control_delay_approx.B, Qy, Rv)
    Kry = lqry(roll_control_delay_approx, Qv, Rv); % continous lqr controller based on output
    Kryd = lqrd(roll_control_delay_approx.A, roll_control_delay_approx.B, Qy, Rv, 0.02 );
    %observer_poles = real(eig(Ari - Bri * Krd)) * 5 + imag(eig(Ari - Bri * Krd))*1i;
    %observer_poles =
    %observer_poles(mod(0:length(roll_control_delay_approx.A)-1, numel(observer_poles)) + 1);
    
    % Observer regulator tuned through pole placement
    observer_poles = [-20:-1:-19 - length(roll_control_delay_approx.A)]';
    L = place(roll_control_delay_approx.A', roll_control_delay_approx.C', observer_poles)';
    controller = reg(roll_control_delay_approx, Kry, L);   
    
    % LQG controller
    kalman_roll_control_delay_approx = minreal(pade(ss(Ari,[Bri Bri],Cri,0, 'InputDelay',[0.07 0.07 0.01 0 0 0]), 5));
    [kest, L, P] = kalman(kalman_roll_control_delay_approx, diag([2 2 2/5]), 0.1);
    kalman_controller = lqgreg(kest, Kry);
    
    Mb_bar = [0.1;0.1;0.01];
    Vd = diag(roll_control_delay_approx.B * Mb_bar);
    Vn = 0.01;
    [L, P, E] = lqe(roll_control_delay_approx.A, roll_control_delay_approx.B, roll_control_delay_approx.C, diag(Mb_bar), Vn * eye(size(roll_control_delay_approx.C, 1)));
    lqe_estimator = estim(roll_control_delay_approx, L, 1:13, 1:3);
  
    lqr(roll_control_delay_approx.A', roll_control_delay_approx.C', Vd, Vn)';
    
    %% BEST ESTIMATOR SO FAR
    L_simulink = lqe(roll_control_delay_approx.A, eye(length(roll_control_delay_approx.A)), roll_control_delay_approx.C, 0.001 * eye(28), 0.0001 * eye(13));
    simulink_estimator = estim(roll_control_delay_approx, L_simulink, 1:13, 1:3);
    
    %{
    for i = 1:10
        try
            kalman(kalman_roll_control_delay_approx, diag([i i i/5]), 0.1);
            i
        catch
        end
    end
    %}
    %estim(roll_control_delay_approx, L, [1:13], [1:3])
    %% No velocity
    q = sym('q', [4 1], 'real');
    qi = sym('gi', [3 1]);
    w = sym('w', [3 1], 'real');
    Mb = sym('Mb', [3 1]);
    I = sym('I', [3 3]);     

    Thrust = sym('Thrust', 'real');
        
    qdot = .5 * quatmultiply(q.', [0;w].').';
    
    derivative_vector = [I \ (Mb - cross(w,I*w)); qdot(2:4); q(2:4)];
    
    Ari = jacobian(derivative_vector, [w; q(2:4); qi]);
    Bri = jacobian(derivative_vector, Mb);
    Cri = eye(size(Ari,1));
    Dri = zeros(size(Bri));
    
    %Inertia = diag([0.0826975856, 0.0826975856, 2.4778e-04]);
    Inertia = diag([0.0829427911790996, 0.0829427911790996, 0.000246795169917015]);
    omega = [0;0;0];

    Ari = subs(Ari, I, Inertia);
    Ari = subs(Ari, w, omega);
    Ari = subs(Ari, q, [1; 0; 0; 0]);
    Ari = subs(Ari, Mb, [0;0;0]);
    
    Bri = subs(Bri, I, Inertia);
    Bri = subs(Bri, q, [1; 0; 0; 0]);
    
    Ari = double(Ari);
    Bri = double(Bri);
    
    Qv = diag([1,1,2, 2,2,0.1, 1,1,0.01]);
    Rv = diag([0.3 0.3 1000]);
    Rd = diag([0.3 0.3 0.6]);

    no_delay = ss(Ari, Bri, Cri, Dri);
    roll_control_delay = ss(Ari, Bri, Cri, Dri, 'InputDelay',[0.07 0.07 0.01]);
    roll_control_delay_approx = minreal(pade(roll_control_delay, 7)); 
    %[a,b,c,d,e] = ctrbf(roll_control_delay_approx.A, roll_control_delay_approx.B, roll_control_delay_approx.C); [sum(e), length(roll_control_delay_approx.A)]
    %[X,K,clp,INFO] = icare(no_delay.A,no_delay.B,Qv,R)
    
    K = lqr(no_delay, Qv, Rv);
    K = lqrd(no_delay.A, no_delay.B, Qv, Rv, 0.02);
    
    L = lqe(roll_control_delay_approx.A, eye(length(roll_control_delay_approx.A)), roll_control_delay_approx.C, 1 * eye(30), eps * eye(9));
            
    K = lqr(roll_control_delay_approx, roll_control_delay_approx.C' * Qv * roll_control_delay_approx.C, Rd);
    K = lqrd(roll_control_delay_approx.A, roll_control_delay_approx.B, roll_control_delay_approx.C' * Qv * roll_control_delay_approx.C, Rd, 0.02);
    sys = reg(roll_control_delay_approx, K, L);
    %{
    [A,B,C,D] = ssdata(no_delay);
                
    nx = size(A, 1);
    ni = 3;
    nu = 3;
    nz = nx+ni;

    C = [eye(ni), zeros(ni, nx - ni)];
    Ai = [A zeros(nx,ni); -C zeros(ni,ni)];
    Bi = [B; -zeros(ni, nu)];
    Ci = [eye(length(Ari)), zeros(length(Ari), length(Ai)-length(Ari));zeros(ni,length(Ai)-ni), eye(ni)];
    Q = Ci' * Qv * Ci;

    K = [];
    i = 0.8;

    while isempty(K) && i > 0
        try
            R = diag([1 1 i]);
            
            [Q,R,N] = ltipack.checkQRS(nx+ni,nu,Q,R,[],{'Q','R','N'});

            % Factor [Q N;N' R] and use square-root formulation when possible
            [F,G,INDEF] = ltipack.factorQRS(Q,R,N);
            if INDEF
               % Proceed with original Q,R,N when [Q N;N' R] is numerically indefinite
               warning(message('Control:design:MustBePositiveDefinite','[Q N;N'' R]','lqi'))
               if Ts==0
                  [X,K,clp,INFO] = icare(Ai,Bi,Q,R,N,E);
               else
                  [X,K,clp,INFO] = idare(Ai,Bi,Q,R,N,E);
               end
            else
                BB = [Bi zeros(nz,nz+nu)];
                QQ = zeros(nz);
                NN = [zeros(nz,nu) F];
                RR = [zeros(nu) G;G' -eye(nz+nu)];
                [X,K,clp,INFO] = icare(Ai,BB,QQ,RR,NN);
                K = K(1:min(nu,end),:);
            end

            switch INFO.Report
               case 1
                  K = [];
                  error('Poor solution accuracy')
               case 2
                  % S and K are not finite
                  error(message('Control:design:lqr1'))
               case 3
                  % Could not compute stabilizing S
                  error(message('Control:design:lqr2'))
            end
        catch
            i = i - 0.001;
        end
    end
    %}
    %% Time invariant MPC controller
    q = sym('q', [4 1], 'real');
    qi = sym('gi', [3 1]);
    w = sym('w', [3 1], 'real');
    Mb = sym('Mb', [3 1]);
    I = sym('I', [3 3]);     
    Ve = sym('Ve', [2 1]);
    Xe = sym('Xe', [2 1]);
    Thrust = sym('Thrust', 'real');
    
    syms q2dcm(qin);
    qin = sym('qin', [1 4], 'real');
    
    q2dcm(qin) = [qin(:,1).^2 + qin(:,2).^2 - qin(:,3).^2 - qin(:,4).^2 , 2.*(qin(:,2).*qin(:,3) + qin(:,1).*qin(:,4)), 2.*(qin(:,2).*qin(:,4) - qin(:,1).*qin(:,3)) ; 2.*(qin(:,2).*qin(:,3) - qin(:,1).*qin(:,4)), qin(:,1).^2 - qin(:,2).^2 + qin(:,3).^2 - qin(:,4).^2, 2.*(qin(:,3).*qin(:,4) + qin(:,1).*qin(:,2)) ; 2.*(qin(:,2).*qin(:,4) + qin(:,1).*qin(:,3)), 2.*(qin(:,3).*qin(:,4) - qin(:,1).*qin(:,2)), qin(:,1).^2 - qin(:,2).^2 - qin(:,3).^2 + qin(:,4).^2];
    
    Inertia = diag([0.0826975856, 0.0826975856, 2.4778e-04]);
    omega = [0;0;0];
    
    Aexyz =  q2dcm(q(1), -q(2), -q(3), -q(4)) * ([0;0;Thrust] + [Mb(1:2);0] / 0.401242753755115); % has to be inverse quat
    qdot = .5 * quatmultiply(q.', [0;w].').';
    
    derivative_vector = [I \ (Mb - cross(w,I*w)); qdot(2:4) ; q(2:4); Aexyz(1:2); Ve];
    Ari = jacobian(derivative_vector, [w; q(2:4); qi; Ve; Xe]);
    Bri = jacobian(derivative_vector, Mb);
    Cri = eye(size(Ari,1));
    Dri = zeros(size(Bri));

    Ari = subs(Ari, I, Inertia);
    Ari = subs(Ari, w, omega);
    Ari = subs(Ari, q, [1; 0; 0; 0]);
    Ari = subs(Ari, Mb, [0;0;0]);
    Ari = subs(Ari, Thrust, 10.6);
    
    Bri = subs(Bri, I, Inertia);
    Bri = subs(Bri, q, [1; 0; 0; 0]);
    
    Ari = double(Ari);
    Bri = double(Bri);
    
    Ts = 0.02;
    time_invariant_sys = ss(Ari, Bri, Cri, Dri, Ts);
    time_invariant_delay_sys = minreal(ss(Ari, Bri, Cri, Dri, 'InputDelay',[0.07 0.07 0.01]));
    discretized = absorbDelay(c2d(minreal(ss(Ari, Bri, Cri, Dri, 'InputDelay',[0.07 0.07 0.01])), 0.02));
    upside_down = mpcstate(mpc1, discretized.C \ [ 0 0 0 1 0 0 0 0 0 0 0 0 0]');
    time_invariant_delay_discrete_sys = minreal(absorbDelay(ss(Ari, Bri, Cri, Dri, 0.02, 'InputDelay',[3 3 1])));
    p = 3;
    m= 1;
    mycobj = mps(sys, Ts, p, m);
    mycobj.MV = struct('Min', -3,'Max',3);
    mpcobj.Weights = struct('MV',0,'MVRate',0.01,'Output',1);
    %% Time varying MPC controller
    q = sym('q', [4 1], 'real');
    qi = sym('gi', [3 1]);
    w = sym('w', [3 1], 'real');
    Mb = sym('Mb', [3 1]);
    I = sym('I', [3 3]);     
    Ve = sym('Ve', [2 1]);
    Xe = sym('Xe', [2 1]);
    Thrust = sym('Thrust', 'real');
    
    syms q2dcm(qin);
    qin = sym('qin', [1 4], 'real');
    
    q2dcm(qin) = [qin(:,1).^2 + qin(:,2).^2 - qin(:,3).^2 - qin(:,4).^2 , 2.*(qin(:,2).*qin(:,3) + qin(:,1).*qin(:,4)), 2.*(qin(:,2).*qin(:,4) - qin(:,1).*qin(:,3)) ; 2.*(qin(:,2).*qin(:,3) - qin(:,1).*qin(:,4)), qin(:,1).^2 - qin(:,2).^2 + qin(:,3).^2 - qin(:,4).^2, 2.*(qin(:,3).*qin(:,4) + qin(:,1).*qin(:,2)) ; 2.*(qin(:,2).*qin(:,4) + qin(:,1).*qin(:,3)), 2.*(qin(:,3).*qin(:,4) - qin(:,1).*qin(:,2)), qin(:,1).^2 - qin(:,2).^2 - qin(:,3).^2 + qin(:,4).^2];
    
    Inertia = diag([0.0826975856, 0.0826975856, 2.4778e-04]);
    omega = [0;0;0];
    
    Aexyz =  q2dcm(q(1), -q(2), -q(3), -q(4)) * ([0;0;Thrust] + [Mb(1:2);0] / 0.401242753755115); % has to be inverse quat
    qdot = .5 * quatmultiply(q.', [0;w].').';
    
    derivative_vector = [I \ (Mb - cross(w,I*w)); qdot(2:4) ; q(2:4); Aexyz(1:2); Ve];
    Ari = jacobian(derivative_vector, [w; q(2:4); qi; Ve; Xe]);
    Bri = jacobian(derivative_vector, Mb);
    Cri = eye(size(Ari,1));
    Dri = zeros(size(Bri));

    Ari = subs(Ari, I, Inertia);
    Ari = subs(Ari, w, omega);
    Ari = subs(Ari, q, [1; 0; 0; 0]);
    Ari = subs(Ari, Mb, [0;0;0]);
    
    Qv = diag([1,1,1, 2,2,2, 1,1,1, 2,2, 1,1]);
    Rv = diag([1 1 0.5]);
    Krv = lqr(A(:,:,1), Bri, Qv, Rv);
    
    %% Time varying LQY
    Thrust_space_lqy = 0.5:2:28.5;
    num_lqr = length(Thrust_space_lqy);
    qz_space = -0.5:0.5:0.5;
    num_qz = length(qz_space);
    qx_space = -1:0.1:1;
    num_qx = length(qx_space);
    qy_space = -1:0.1:1;
    num_qy = length(qy_space);
    ss_lqr = sampleBlock(tunable_full, 'Thrust', 10.6 , 'Ixx', 0.0826975856, 'Izz', 2.4778e-04, 'q', [1 0 0 0], 'w', [0 0 0], 'Mb', [0 0 0]);
   
    Qy = double(tunable_full.C)' * diag([1,1,1, 2,2,2, 1,1,1, 2,2, 1,1]) * double(tunable_full.C);
    %% Time varying LQR
    num_lqr = 28;
    Thrust_space_lqr = linspace(0.5,28.5,num_lqr);
    ss_lqr = sampleBlock(tunable, 'T', Thrust_space_lqr);
    
    Qv = diag([1,1,1, 2,2,2, 1,1,1, 2,2, 1,1]);
    %% LQI
    q = sym('q', [4 1], 'real');
    w = sym('w', [3 1], 'real');
    Mb = sym('Mb', [3 1]);
    I = sym('I', [3 3]);     
    Ve = sym('Ve', [2 1]);
    Thrust = sym('Thrust', 'real');
 
    syms q2dcm(qin);
    qin = sym('qin', [1 4], 'real');
    
    q2dcm(qin) = [qin(:,1).^2 + qin(:,2).^2 - qin(:,3).^2 - qin(:,4).^2 , 2.*(qin(:,2).*qin(:,3) + qin(:,1).*qin(:,4)), 2.*(qin(:,2).*qin(:,4) - qin(:,1).*qin(:,3)) ; 2.*(qin(:,2).*qin(:,3) - qin(:,1).*qin(:,4)), qin(:,1).^2 - qin(:,2).^2 + qin(:,3).^2 - qin(:,4).^2, 2.*(qin(:,3).*qin(:,4) + qin(:,1).*qin(:,2)) ; 2.*(qin(:,2).*qin(:,4) + qin(:,1).*qin(:,3)), 2.*(qin(:,3).*qin(:,4) - qin(:,1).*qin(:,2)), qin(:,1).^2 - qin(:,2).^2 - qin(:,3).^2 + qin(:,4).^2];
    
    Aexyz =  q2dcm(q(1), -q(2), -q(3), -q(4)) * ([0;0;Thrust] + [Mb(1:2);0] / 0.401242753755115); % has to be inverse quat
    qdot = .5 * quatmultiply(q.', [0;w].').';
    
    derivative_vector = [I \ (Mb - cross(w,I*w)); qdot(2:4); Aexyz(1:2)];
    
    Ari = jacobian(derivative_vector, [w; q(2:4); Ve]);
    Bri = jacobian(derivative_vector, Mb);
    Cri = eye(size(Ari,1));
    Dri = zeros(size(Bri));
    
    Inertia = diag([0.0826975856, 0.0826975856, 2.4778e-04]);
    omega = [0;0;0];

    Ari = subs(Ari, I, Inertia);
    Ari = subs(Ari, w, omega);
    Ari = subs(Ari, Mb, [0;0;0]);
    Ari = subs(Ari, q, [1 0 0 0]');
    Ari = subs(Ari, Thrust, 10.6);
    
    Bri = subs(Bri, I, Inertia);
    Bri = subs(Bri, q, [1 0 0 0]');
    
    Ari = double(Ari);
    Bri = double(Bri);
    
    roll_control_lqi = ss(Ari, Bri, Cri, Dri, 'InputDelay', [0.07 0.07 0.01]);
    roll_control_delay_lqi_pade = minreal(pade(roll_control_lqi, 5));
    %[a,b,c,d,e] = ctrbf(roll_control_delay_lqi_pade.A, roll_control_delay_lqi_pade.B, roll_control_delay_lqi_pade.C); [sum(e), length(roll_control_delay_lqi_pade.A)]
    Qy = diag([1,1,1, 2,2,2, 2,2, 2,2,1, 1,1]);
    
    [A,B,C,D] = ssdata(roll_control_delay_lqi_pade);
                
    nx = size(A, 1);
    ni = 5;
    nu = 3;
    nz = nx+ni;

    C = [eye(ni), zeros(ni, nx - ni)];
    Ai = [A zeros(nx,ni); -C zeros(ni,ni)];
    Bi = [B; -zeros(ni, nu)];
    Ci = [eye(8), zeros(8, length(Ai)-8);zeros(5,length(Ai)-5), eye(5)];
    Q = Ci' * Qy * Ci;

    sys = minreal(ss(Ai, Bi, Ci, 0));
    
    Ai = sys.A;
    Bi = sys.B;
    Ci = sys.C;
    Q = Ci' * Qy * Ci;
    nx = size(Ai, 1) - 5;
    nz = size(Ai, 1);
    %{
    BB = [Bi zeros(nz,nz+nu)];
    [a,b,c,d,e] = ctrbf(Ai, BB, Ci); 
    [sum(e), length(Ai)]
    %}
    %{
    Ci = [eye(8), zeros(8, 20);zeros(5,23), eye(5)];
    BB = [Bi zeros(nz,nz+nu)];
    sys = minreal(ss(Ai, BB, Ci, 0));
    
    Ai = sys.A;
    Bi = sys.B;
    Ci = sys.C;

    [a,b,c,d,e] = ctrbf(Ai, Bi, Ci); 
    [sum(e), length(Ai)]
    %}
    K = [];
    i = 1;

    while isempty(K) && i > 0
        try
            R = diag([0.3 0.3 i]);
            
            [Q,R,N] = ltipack.checkQRS(nx+ni,nu,Q,R,[],{'Q','R','N'});

            % Factor [Q N;N' R] and use square-root formulation when possible
            [F,G,INDEF] = ltipack.factorQRS(Q,R,N);
            if INDEF
               % Proceed with original Q,R,N when [Q N;N' R] is numerically indefinite
               warning(message('Control:design:MustBePositiveDefinite','[Q N;N'' R]','lqi'))
               if Ts==0
                  [X,K,clp,INFO] = icare(Ai,Bi,Q,R,N,E);
               else
                  [X,K,clp,INFO] = idare(Ai,Bi,Q,R,N,E);
               end
            else
                BB = [Bi zeros(nz,nz+nu)];
                QQ = zeros(nz);
                NN = [zeros(nz,nu) F];
                RR = [zeros(nu) G;G' -eye(nz+nu)];
                [X,K,clp,INFO] = icare(Ai,BB,QQ,RR,NN);
                K = K(1:min(nu,end),:);
            end

            switch INFO.Report
               case 1
                  K = [];
                  error('Poor solution accuracy')
               case 2
                  % S and K are not finite
                  error(message('Control:design:lqr1'))
               case 3
                  % Could not compute stabilizing S
                  error(message('Control:design:lqr2'))
            end
        catch
            i = i - 0.001;
        end
    end
    
    l = 0.001;
    clear L
    while ~exist('L', 'var')
        try
             L = lqe(roll_control_delay_approx.A, eye(length(roll_control_delay_approx.A)), roll_control_delay_approx.C, l * eye(28), 0.0001 * eye(13));
        catch
            i = i - 0.00001;
        end
    end

    sys = reg(roll_control_delay_approx, Kry, L);
    %% tunable state space
    q = realp('q', [1 0 0 0]);
    qi = realp('gi', [0 0 0]);
    w = realp('w', [0 0 0]);
    Mb = realp('Mb', [0 0 0]);
    Ixx = realp('Ixx', 0.0826975856);
    Izz = realp('Izz', 2.4778e-04);
    
    I = [Ixx 0 0; 0 Ixx 0; 0 0 Izz];
    Ve = realp('Ve', [0 0]);
    Xe = realp('Xe', [0 0]);
    Thrust = realp('Thrust', 10.6);
    
    A = cat(2, [(I(1,3)*I(2,2)^2*w(2) - I(1,2)*I(2,3)^2*w(3) + I(1,3)*I(3,2)^2*w(2) - I(1,2)*I(3,3)^2*w(3) + I(1,1)*I(1,2)*I(2,3)*w(2) - I(1,1)*I(1,3)*I(2,2)*w(2) - 2*I(1,2)*I(2,1)*I(2,3)*w(1) + 2*I(1,3)*I(2,1)*I(2,2)*w(1) + I(1,1)*I(1,2)*I(3,3)*w(3) - I(1,1)*I(1,3)*I(3,2)*w(3) - I(1,2)*I(2,2)*I(2,3)*w(2) + I(1,3)*I(2,2)*I(2,3)*w(3) - 2*I(1,2)*I(3,1)*I(3,3)*w(1) + 2*I(1,3)*I(3,1)*I(3,2)*w(1) - I(1,2)*I(3,2)*I(3,3)*w(2) + I(2,1)*I(2,2)*I(3,3)*w(3) - I(2,1)*I(2,3)*I(3,2)*w(3) + I(1,3)*I(3,2)*I(3,3)*w(3) - I(2,2)*I(3,1)*I(3,3)*w(2) + I(2,3)*I(3,1)*I(3,2)*w(2))/(I(1,1)*I(2,2)*I(3,3) - I(1,1)*I(2,3)*I(3,2) - I(1,2)*I(2,1)*I(3,3) + I(1,2)*I(2,3)*I(3,1) + I(1,3)*I(2,1)*I(3,2) - I(1,3)*I(2,2)*I(3,1)), (w(3)*I(1,2)^2*I(3,3) + 2*I(2,3)*w(2)*I(1,2)^2 - 2*w(2)*I(1,2)*I(1,3)*I(2,2) - w(3)*I(1,2)*I(1,3)*I(3,2) + I(2,3)*w(3)*I(1,2)*I(1,3) - I(2,3)*w(1)*I(1,2)*I(2,2) - w(1)*I(1,2)*I(3,2)*I(3,3) + I(1,1)*I(2,3)*w(1)*I(1,2) - w(3)*I(1,3)^2*I(2,2) + w(1)*I(1,3)*I(2,2)^2 - I(1,1)*w(1)*I(1,3)*I(2,2) + w(1)*I(1,3)*I(3,2)^2 + w(3)*I(2,2)^2*I(3,3) - 2*w(2)*I(2,2)*I(3,2)*I(3,3) - I(2,3)*w(3)*I(2,2)*I(3,2) - w(3)*I(2,2)*I(3,3)^2 - I(3,1)*w(1)*I(2,2)*I(3,3) + 2*I(2,3)*w(2)*I(3,2)^2 + I(2,3)*w(3)*I(3,2)*I(3,3) + I(2,3)*I(3,1)*w(1)*I(3,2))/(I(1,1)*I(2,2)*I(3,3) - I(1,1)*I(2,3)*I(3,2) - I(1,2)*I(2,1)*I(3,3) + I(1,2)*I(2,3)*I(3,1) + I(1,3)*I(2,1)*I(3,2) - I(1,3)*I(2,2)*I(3,1)), -(- w(2)*I(1,2)^2*I(3,3) - w(2)*I(1,2)*I(1,3)*I(2,3) - 2*w(3)*I(1,2)*I(1,3)*I(3,3) + I(3,2)*w(2)*I(1,2)*I(1,3) + w(1)*I(1,2)*I(2,3)^2 + w(1)*I(1,2)*I(3,3)^2 - I(1,1)*w(1)*I(1,2)*I(3,3) + w(2)*I(1,3)^2*I(2,2) + 2*I(3,2)*w(3)*I(1,3)^2 - w(1)*I(1,3)*I(2,2)*I(2,3) - I(3,2)*w(1)*I(1,3)*I(3,3) + I(1,1)*I(3,2)*w(1)*I(1,3) - w(2)*I(2,2)^2*I(3,3) - 2*w(3)*I(2,2)*I(2,3)*I(3,3) + I(3,2)*w(2)*I(2,2)*I(2,3) + w(2)*I(2,2)*I(3,3)^2 - I(2,1)*w(1)*I(2,2)*I(3,3) + 2*I(3,2)*w(3)*I(2,3)^2 - I(3,2)*w(2)*I(2,3)*I(3,3) + I(2,1)*I(3,2)*w(1)*I(2,3))/(I(1,1)*I(2,2)*I(3,3) - I(1,1)*I(2,3)*I(3,2) - I(1,2)*I(2,1)*I(3,3) + I(1,2)*I(2,3)*I(3,1) + I(1,3)*I(2,1)*I(3,2) - I(1,3)*I(2,2)*I(3,1)), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
     -(w(2)*I(1,1)^2*I(2,3) + w(3)*I(1,1)^2*I(3,3) - 2*w(1)*I(1,1)*I(2,1)*I(2,3) - I(1,3)*w(2)*I(1,1)*I(2,1) - w(3)*I(1,1)*I(2,3)^2 - I(2,2)*w(2)*I(1,1)*I(2,3) - 2*w(1)*I(1,1)*I(3,1)*I(3,3) - I(1,3)*w(3)*I(1,1)*I(3,1) - w(3)*I(1,1)*I(3,3)^2 - I(3,2)*w(2)*I(1,1)*I(3,3) + w(3)*I(2,1)^2*I(3,3) + 2*I(1,3)*w(1)*I(2,1)^2 - w(3)*I(2,1)*I(2,3)*I(3,1) + I(1,3)*w(3)*I(2,1)*I(2,3) - w(2)*I(2,1)*I(3,1)*I(3,3) + I(1,3)*I(2,2)*w(2)*I(2,1) + w(2)*I(2,3)*I(3,1)^2 + 2*I(1,3)*w(1)*I(3,1)^2 + I(1,3)*w(3)*I(3,1)*I(3,3) + I(1,3)*I(3,2)*w(2)*I(3,1))/(I(1,1)*I(2,2)*I(3,3) - I(1,1)*I(2,3)*I(3,2) - I(1,2)*I(2,1)*I(3,3) + I(1,2)*I(2,3)*I(3,1) + I(1,3)*I(2,1)*I(3,2) - I(1,3)*I(2,2)*I(3,1)), -(I(1,1)^2*I(2,3)*w(1) - I(1,3)^2*I(2,1)*w(3) + I(2,3)*I(3,1)^2*w(1) - I(2,1)*I(3,3)^2*w(3) - I(1,1)*I(1,3)*I(2,1)*w(1) + 2*I(1,1)*I(1,2)*I(2,3)*w(2) - 2*I(1,2)*I(1,3)*I(2,1)*w(2) + I(1,1)*I(1,3)*I(2,3)*w(3) - I(1,1)*I(2,2)*I(2,3)*w(1) + I(1,3)*I(2,1)*I(2,2)*w(1) + I(1,1)*I(1,2)*I(3,3)*w(3) - I(1,2)*I(1,3)*I(3,1)*w(3) - I(1,1)*I(3,2)*I(3,3)*w(1) + I(1,3)*I(3,1)*I(3,2)*w(1) + I(2,1)*I(2,2)*I(3,3)*w(3) - I(2,2)*I(2,3)*I(3,1)*w(3) - I(2,1)*I(3,1)*I(3,3)*w(1) - 2*I(2,1)*I(3,2)*I(3,3)*w(2) + 2*I(2,3)*I(3,1)*I(3,2)*w(2) + I(2,3)*I(3,1)*I(3,3)*w(3))/(I(1,1)*I(2,2)*I(3,3) - I(1,1)*I(2,3)*I(3,2) - I(1,2)*I(2,1)*I(3,3) + I(1,2)*I(2,3)*I(3,1) + I(1,3)*I(2,1)*I(3,2) - I(1,3)*I(2,2)*I(3,1)), (- w(1)*I(1,1)^2*I(3,3) - w(2)*I(1,1)*I(1,3)*I(2,3) - 2*w(3)*I(1,1)*I(1,3)*I(3,3) + I(3,1)*w(1)*I(1,1)*I(1,3) + w(1)*I(1,1)*I(2,3)^2 + w(1)*I(1,1)*I(3,3)^2 - I(1,2)*w(2)*I(1,1)*I(3,3) + w(2)*I(1,3)^2*I(2,1) + 2*I(3,1)*w(3)*I(1,3)^2 - w(1)*I(1,3)*I(2,1)*I(2,3) - I(3,1)*w(1)*I(1,3)*I(3,3) + I(1,2)*I(3,1)*w(2)*I(1,3) - w(1)*I(2,1)^2*I(3,3) - 2*w(3)*I(2,1)*I(2,3)*I(3,3) + I(3,1)*w(1)*I(2,1)*I(2,3) + w(2)*I(2,1)*I(3,3)^2 - I(2,2)*w(2)*I(2,1)*I(3,3) + 2*I(3,1)*w(3)*I(2,3)^2 - I(3,1)*w(2)*I(2,3)*I(3,3) + I(2,2)*I(3,1)*w(2)*I(2,3))/(I(1,1)*I(2,2)*I(3,3) - I(1,1)*I(2,3)*I(3,2) - I(1,2)*I(2,1)*I(3,3) + I(1,2)*I(2,3)*I(3,1) + I(1,3)*I(2,1)*I(3,2) - I(1,3)*I(2,2)*I(3,1)), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
     (w(2)*I(1,1)^2*I(2,2) + w(3)*I(1,1)^2*I(3,2) - 2*w(1)*I(1,1)*I(2,1)*I(2,2) - I(1,2)*w(2)*I(1,1)*I(2,1) - w(2)*I(1,1)*I(2,2)^2 - I(2,3)*w(3)*I(1,1)*I(2,2) - 2*w(1)*I(1,1)*I(3,1)*I(3,2) - I(1,2)*w(3)*I(1,1)*I(3,1) - w(2)*I(1,1)*I(3,2)^2 - I(3,3)*w(3)*I(1,1)*I(3,2) + w(3)*I(2,1)^2*I(3,2) + 2*I(1,2)*w(1)*I(2,1)^2 - w(3)*I(2,1)*I(2,2)*I(3,1) + I(1,2)*w(2)*I(2,1)*I(2,2) - w(2)*I(2,1)*I(3,1)*I(3,2) + I(1,2)*I(2,3)*w(3)*I(2,1) + w(2)*I(2,2)*I(3,1)^2 + 2*I(1,2)*w(1)*I(3,1)^2 + I(1,2)*w(2)*I(3,1)*I(3,2) + I(1,2)*I(3,3)*w(3)*I(3,1))/(I(1,1)*I(2,2)*I(3,3) - I(1,1)*I(2,3)*I(3,2) - I(1,2)*I(2,1)*I(3,3) + I(1,2)*I(2,3)*I(3,1) + I(1,3)*I(2,1)*I(3,2) - I(1,3)*I(2,2)*I(3,1)), -(- w(1)*I(1,1)^2*I(2,2) - 2*w(2)*I(1,1)*I(1,2)*I(2,2) - w(3)*I(1,1)*I(1,2)*I(3,2) + I(2,1)*w(1)*I(1,1)*I(1,2) + w(1)*I(1,1)*I(2,2)^2 - I(1,3)*w(3)*I(1,1)*I(2,2) + w(1)*I(1,1)*I(3,2)^2 + w(3)*I(1,2)^2*I(3,1) + 2*I(2,1)*w(2)*I(1,2)^2 - I(2,1)*w(1)*I(1,2)*I(2,2) - w(1)*I(1,2)*I(3,1)*I(3,2) + I(1,3)*I(2,1)*w(3)*I(1,2) + w(3)*I(2,2)^2*I(3,1) - w(1)*I(2,2)*I(3,1)^2 - 2*w(2)*I(2,2)*I(3,1)*I(3,2) - I(3,3)*w(3)*I(2,2)*I(3,1) - I(2,1)*w(3)*I(2,2)*I(3,2) + I(2,1)*w(1)*I(3,1)*I(3,2) + 2*I(2,1)*w(2)*I(3,2)^2 + I(2,1)*I(3,3)*w(3)*I(3,2))/(I(1,1)*I(2,2)*I(3,3) - I(1,1)*I(2,3)*I(3,2) - I(1,2)*I(2,1)*I(3,3) + I(1,2)*I(2,3)*I(3,1) + I(1,3)*I(2,1)*I(3,2) - I(1,3)*I(2,2)*I(3,1)), (I(1,1)^2*I(3,2)*w(1) - I(1,2)^2*I(3,1)*w(2) + I(2,1)^2*I(3,2)*w(1) - I(2,2)^2*I(3,1)*w(2) + I(1,1)*I(1,3)*I(2,2)*w(2) - I(1,2)*I(1,3)*I(2,1)*w(2) - I(1,1)*I(1,2)*I(3,1)*w(1) + I(1,1)*I(1,2)*I(3,2)*w(2) - I(1,1)*I(2,2)*I(2,3)*w(1) + I(1,2)*I(2,1)*I(2,3)*w(1) + 2*I(1,1)*I(1,3)*I(3,2)*w(3) - 2*I(1,2)*I(1,3)*I(3,1)*w(3) - I(2,1)*I(2,2)*I(3,1)*w(1) - I(1,1)*I(3,2)*I(3,3)*w(1) + I(1,2)*I(3,1)*I(3,3)*w(1) + I(2,1)*I(2,2)*I(3,2)*w(2) + 2*I(2,1)*I(2,3)*I(3,2)*w(3) - 2*I(2,2)*I(2,3)*I(3,1)*w(3) - I(2,1)*I(3,2)*I(3,3)*w(2) + I(2,2)*I(3,1)*I(3,3)*w(2))/(I(1,1)*I(2,2)*I(3,3) - I(1,1)*I(2,3)*I(3,2) - I(1,2)*I(2,1)*I(3,3) + I(1,2)*I(2,3)*I(3,1) + I(1,3)*I(2,1)*I(3,2) - I(1,3)*I(2,2)*I(3,1)), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
     q(1)/2, -q(4)/2, q(3)/2, 0, w(3)/2, -w(2)/2, 0, 0, 0, 0, 0, 0, 0;
     q(4)/2, q(1)/2, -q(2)/2, -w(3)/2, 0, w(1)/2, 0, 0, 0, 0, 0, 0, 0;
     -q(3)/2, q(2)/2, q(1)/2, w(2)/2, -w(1)/2, 0, 0, 0, 0, 0, 0, 0, 0;
     0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0;
     0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0;
     0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0;
     0, 0, 0, (18014398509481984*Mb(1)*q(2))/3614073432593295 + (18014398509481984*Mb(2)*q(3))/3614073432593295 + 2*Thrust*q(4), (18014398509481984*Mb(2)*q(2))/3614073432593295 - (18014398509481984*Mb(1)*q(3))/3614073432593295 + 2*Thrust*q(1), 2*Thrust*q(2) - (18014398509481984*Mb(1)*q(4))/3614073432593295 - (18014398509481984*Mb(2)*q(1))/3614073432593295, 0, 0, 0, 0, 0, 0, 0;
     0, 0, 0, (18014398509481984*Mb(1)*q(3))/3614073432593295 - (18014398509481984*Mb(2)*q(2))/3614073432593295 - 2*Thrust*q(1), (18014398509481984*Mb(1)*q(2))/3614073432593295 + (18014398509481984*Mb(2)*q(3))/3614073432593295 + 2*Thrust*q(4), (18014398509481984*Mb(1)*q(1))/3614073432593295 - (18014398509481984*Mb(2)*q(4))/3614073432593295 + 2*Thrust*q(3), 0, 0, 0, 0, 0, 0, 0;
     0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0;
     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;
     0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,0,0], [323.899592610528,0,112.967212824542,0,25.3284666827986,0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,323.899592610528,0,112.967212824542,0,25.3284666827986,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,189179.917668900,0,101032.976049980,0,17343.4597864117;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     66.7568798066371,0,23.2829519414889,0,5.22028877036333,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,66.7568798066371,0,23.2829519414889,0,5.22028877036333,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     -428.571428571429,-334.821428571429,-149.473852040816,-75.0705730115707,-33.5136486658798,0,0,0,0,0,0,0,0,0,0;
     256,0,0,0,0,0,0,0,0,0,0,0,0,0,0;0,256,0,0,0,0,0,0,0,0,0,0,0,0,0;0,0,128,0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,64,0,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,-428.571428571429,-334.821428571429,-149.473852040816,-75.0705730115707,-33.5136486658798,0,0,0,0,0;
     0,0,0,0,0,256,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,256,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,128,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,64,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,-3000.00000000000,-2050.78125000000,-1602.17285156250,-704.079866409302,-275.031197816133;0,0,0,0,0,0,0,0,0,0,2048,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,1024,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,0,1024,0,0;
     0,0,0,0,0,0,0,0,0,0,0,0,0,512,0]);

    B = [(I(2,2)*I(3,3) - I(2,3)*I(3,2))/(I(1,1)*I(2,2)*I(3,3) - I(1,1)*I(2,3)*I(3,2) - I(1,2)*I(2,1)*I(3,3) + I(1,2)*I(2,3)*I(3,1) + I(1,3)*I(2,1)*I(3,2) - I(1,3)*I(2,2)*I(3,1)), -(I(1,2)*I(3,3) - I(1,3)*I(3,2))/(I(1,1)*I(2,2)*I(3,3) - I(1,1)*I(2,3)*I(3,2) - I(1,2)*I(2,1)*I(3,3) + I(1,2)*I(2,3)*I(3,1) + I(1,3)*I(2,1)*I(3,2) - I(1,3)*I(2,2)*I(3,1)), (I(1,2)*I(2,3) - I(1,3)*I(2,2))/(I(1,1)*I(2,2)*I(3,3) - I(1,1)*I(2,3)*I(3,2) - I(1,2)*I(2,1)*I(3,3) + I(1,2)*I(2,3)*I(3,1) + I(1,3)*I(2,1)*I(3,2) - I(1,3)*I(2,2)*I(3,1));
     -(I(2,1)*I(3,3) - I(2,3)*I(3,1))/(I(1,1)*I(2,2)*I(3,3) - I(1,1)*I(2,3)*I(3,2) - I(1,2)*I(2,1)*I(3,3) + I(1,2)*I(2,3)*I(3,1) + I(1,3)*I(2,1)*I(3,2) - I(1,3)*I(2,2)*I(3,1)), (I(1,1)*I(3,3) - I(1,3)*I(3,1))/(I(1,1)*I(2,2)*I(3,3) - I(1,1)*I(2,3)*I(3,2) - I(1,2)*I(2,1)*I(3,3) + I(1,2)*I(2,3)*I(3,1) + I(1,3)*I(2,1)*I(3,2) - I(1,3)*I(2,2)*I(3,1)), -(I(1,1)*I(2,3) - I(1,3)*I(2,1))/(I(1,1)*I(2,2)*I(3,3) - I(1,1)*I(2,3)*I(3,2) - I(1,2)*I(2,1)*I(3,3) + I(1,2)*I(2,3)*I(3,1) + I(1,3)*I(2,1)*I(3,2) - I(1,3)*I(2,2)*I(3,1));
     (I(2,1)*I(3,2) - I(2,2)*I(3,1))/(I(1,1)*I(2,2)*I(3,3) - I(1,1)*I(2,3)*I(3,2) - I(1,2)*I(2,1)*I(3,3) + I(1,2)*I(2,3)*I(3,1) + I(1,3)*I(2,1)*I(3,2) - I(1,3)*I(2,2)*I(3,1)), -(I(1,1)*I(3,2) - I(1,2)*I(3,1))/(I(1,1)*I(2,2)*I(3,3) - I(1,1)*I(2,3)*I(3,2) - I(1,2)*I(2,1)*I(3,3) + I(1,2)*I(2,3)*I(3,1) + I(1,3)*I(2,1)*I(3,2) - I(1,3)*I(2,2)*I(3,1)), (I(1,1)*I(2,2) - I(1,2)*I(2,1))/(I(1,1)*I(2,2)*I(3,3) - I(1,1)*I(2,3)*I(3,2) - I(1,2)*I(2,1)*I(3,3) + I(1,2)*I(2,3)*I(3,1) + I(1,3)*I(2,1)*I(3,2) - I(1,3)*I(2,2)*I(3,1));
     0, 0, 0;
     0, 0, 0;
     0, 0, 0;
     0, 0, 0;
     0, 0, 0;
     0, 0, 0;
     (9007199254740992*q(1)^2)/3614073432593295 + (9007199254740992*q(2)^2)/3614073432593295 - (9007199254740992*q(3)^2)/3614073432593295 - (9007199254740992*q(4)^2)/3614073432593295, (18014398509481984*q(2)*q(3))/3614073432593295 - (18014398509481984*q(1)*q(4))/3614073432593295, 0;
     (18014398509481984*q(1)*q(4))/3614073432593295 + (18014398509481984*q(2)*q(3))/3614073432593295, (9007199254740992*q(1)^2)/3614073432593295 + (9007199254740992*q(3)^2)/3614073432593295 - (9007199254740992*q(2)^2)/3614073432593295 - (9007199254740992*q(4)^2)/3614073432593295, 0;
     0, 0, 0;
     0, 0, 0;
     32,0,0;
     0,0,0;
     0,0,0;
     0,0,0;
     0,0,0;
     0,32,0;
     0,0,0;
     0,0,0;
     0,0,0;
     0,0,0;
     0,0,128;
     0,0,0;
     0,0,0;
     0,0,0;
     0,0,0];

    C = cat(2, eye(13), zeros(13, 15));
    D = zeros(13, 3);
    tunable_full = ss(A, B, C, D);

    T = realp('T', 10.6);
    A_tunable = [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
                 1/2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
                 0, 1/2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
                 0, 0, 1/2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
                 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0;
                 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0;
                 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0;
                 0, 0, 0, 0, 2*T, 0, 0, 0, 0, 0, 0, 0, 0;
                 0, 0, 0, -2*T, 0, 0, 0, 0, 0, 0, 0, 0, 0;
                 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0;
                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0];
    tunable = ss(A_tunable, Bri, Cri, Dri);
    
    Thrust_space = 1:28;
    ss_sample = sampleBlock(tunable, 'T', Thrust_space);
    ss_sample.SamplingGrid = struct('T', Thrust_space);
    domain = struct('T', Thrust_space.');
    shapefcn = polyBasis('canonical',1);
    Kinit = [1.5358    0.3831    0.0000   15.4029    0.4122    0.0000    0.0000   -0.0000    0.0000    0.0000   -1.8589    0.0000   -1.0000;
   -0.3831    1.5358   -0.0000   -0.4122   15.4029   -0.0000   -0.0000   -0.0000   -0.0000    1.8589   -0.0000    1.0000    0.0000;
   -0.0000    0.0000    1.2913   -0.0000    0.0000    3.1625   -0.0000   -0.0000    1.2910    0.0000   -0.0000    0.0000    0.0000];
    K_tunable = tunableSurface('K', Kinit, domain, shapefcn);
    
    open('Gain_Scheduler');
    ST = slTuner('Gain_Scheduler', 'K');
    ST.addPoint({'Disturbance' 'y'});
    ST.setBlockParam('K', K_tunable);
    
    Rrejection = TuningGoal.StepRejection('Disturbance', 'y(1)', 0.1, 2);
    ST = systune(ST,Rrejection)
    writeBlockValue(ST);
    
    %% No quaternion integral, velocity control
    q = sym('q', [4 1], 'real');
    w = sym('w', [3 1], 'real');
    Mb = sym('Mb', [3 1]);
    I = sym('I', [3 3]);     
    Ve = sym('Ve', [2 1]);
    Xe = sym('Xe', [2 1]);
    Thrust = sym('Thrust', 'real');

    syms q2dcm(qin);
    qin = sym('qin', [1 4], 'real');

    q2dcm(qin) = [qin(:,1).^2 + qin(:,2).^2 - qin(:,3).^2 - qin(:,4).^2 , 2.*(qin(:,2).*qin(:,3) + qin(:,1).*qin(:,4)), 2.*(qin(:,2).*qin(:,4) - qin(:,1).*qin(:,3)) ; 2.*(qin(:,2).*qin(:,3) - qin(:,1).*qin(:,4)), qin(:,1).^2 - qin(:,2).^2 + qin(:,3).^2 - qin(:,4).^2, 2.*(qin(:,3).*qin(:,4) + qin(:,1).*qin(:,2)) ; 2.*(qin(:,2).*qin(:,4) + qin(:,1).*qin(:,3)), 2.*(qin(:,3).*qin(:,4) - qin(:,1).*qin(:,2)), qin(:,1).^2 - qin(:,2).^2 - qin(:,3).^2 + qin(:,4).^2];

    Aexyz =  q2dcm(q(1), -q(2), -q(3), -q(4)) * ([0;0;Thrust] + [Mb(1:2);0] / 0.401242753755115); % has to be inverse quat
    qdot = .5 * quatmultiply(q.', [0;w].').';

    derivative_vector = [I \ (Mb - cross(w,I*w)); qdot(2:4); Aexyz(1:2); Ve];

    Ari = jacobian(derivative_vector, [w; q(2:4); Ve; Xe]);
    Bri = jacobian(derivative_vector, Mb);
    Cri = eye(size(Ari,1));
    Dri = zeros(size(Bri));

    Inertia = diag([0.0826975856, 0.0826975856, 2.4778e-04]);
    omega = [0;0;0];

    Ari = subs(Ari, I, Inertia);
    Ari = subs(Ari, w, omega);
    Ari = subs(Ari, q, [1; 0; 0; 0]);
    Ari = subs(Ari, Mb, [0;0;0]);
    Ari = subs(Ari, Thrust, 10.6);

    Bri = subs(Bri, I, Inertia);
    Bri = subs(Bri, q, [1; 0; 0; 0]);

    Ari = double(Ari);
    Bri = double(Bri);
    %Bref = [Bri, [zeros(6,3);-eye(3);zeros(4,3)]];

    Qv = diag([1,1,1, 2,2,2, 1,1, 1,1]);
    Rv = diag([0.3 0.3 0.6]);
    Krv = lqr(Ari, Bri, Qv, diag([0.3 0.3 0.47]) );
    %Krv(abs(Krv)<1e-5)=0;

    Krd = lqrd(Ari, Bri, Qv, Rv, 0.02); % discrete lqr controller for continous system

    roll_control_v = ss(Ari, Bri, Cri, Dri);
    roll_control_discrete = c2d(roll_control_v, 0.02);
    roll_control_delay = ss(Ari, Bri, Cri, Dri, 'InputDelay',[0.07 0.07 0.01]);
    roll_control_delay_approx = minreal(pade(roll_control_delay, 3)); 
    %[a,b,c,d,e] = ctrbf(roll_control_delay_approx.A, roll_control_delay_approx.B, roll_control_delay_approx.C); [sum(e), length(roll_control_delay_approx.A)]

    Qy = roll_control_delay_approx.C' * Qv * roll_control_delay_approx.C;
    % lqr(roll_control_delay_approx, Qy, Rv), same as lqry
    %     icare(roll_control_delay_approx.A, roll_control_delay_approx.B, Qy, Rv)
    Kry = lqry(roll_control_delay_approx, Qv, Rv); % continous lqr controller based on output
    Kryd = lqrd(roll_control_delay_approx.A, roll_control_delay_approx.B, Qy, Rv, 0.02 );
    %observer_poles = real(eig(Ari - Bri * Krd)) * 5 + imag(eig(Ari - Bri * Krd))*1i;
    %observer_poles =
    %observer_poles(mod(0:length(roll_control_delay_approx.A)-1, numel(observer_poles)) + 1);

    % Observer regulator tuned through pole placement
    observer_poles = [-20:-1:-19 - length(roll_control_delay_approx.A)]';
    L = place(roll_control_delay_approx.A', roll_control_delay_approx.C', observer_poles)';
    controller = reg(roll_control_delay_approx, Kry, L);   

    % LQG controllers
    [kest, L, P] = kalman(roll_control_delay_approx, diag([2 2 2/5]), 0.1 * eye(size(roll_control_delay_approx.C, 1)));
    kalman_controller = lqgreg(kest, Kry);

    Mb_bar = [0.1;0.1;0.01];
    Vd = diag(roll_control_delay_approx.B * Mb_bar);
    Vn = 0.01 * eye(size(roll_control_delay_approx.C, 1));
    [L, P, E] = lqe(roll_control_delay_approx.A, roll_control_delay_approx.B, roll_control_delay_approx.C, diag(Mb_bar), Vn);
    lqe_estimator = estim(roll_control_delay_approx, L, 1:size(roll_control_delay_approx.C, 1), 1:3);

    % lqr(roll_control_delay_approx.A', roll_control_delay_approx.C', Vd, Vn)';  