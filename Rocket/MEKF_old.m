function MEKF_old
    dt = 0.001;
    tstart = -5;
    time = 30;
    
    I = diag([0.08 0.08 0.002]);
    % Mb = [0.05; -0.02; 0.00005];
    P = 1e-10 * eye(18);
    w = [0; 0; 0];
    Ve = [0; 0; 0];
    Xe = [0; 0; 0];
    q = [1 0 0 0];
    nominal_state = [0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0].';
    g = 9.08665;

    % error_w = zeros(4, time/dt);
    % store_nominal_state = zeros(19, time/dt);

    gyro_bias_sym = sym('gyro_bias', [3 1]);
    accelerometer_bias_sym = sym('accelerometer_bias', [3 1]);
    accelerometer_true_measurements_sym = sym('accelerometer_true_measurements', [3 1]);
    angular_velocity_sym = sym('angular_velocity', [3 1]);
    Xe_sym = sym('Xe_sym', [3 1]);
    Ve_sym = sym('Ve_sym', [3 1]);
    quat_sym = sym('quat_sym', [1 4], 'real');
    sea_level_pressure_sym = sym('sea_level_pressure_sym', 1);
    temperature = sym('temperature', 1);
    
    delta_Xe = sym('delta_Xe', [3 1]);
    delta_Ve = sym('delta_Ve', [3 1]);
    delta_theta = sym('delta_theta', [3 1]);
    delta_ab = sym('delta_ab', [3 1]);
    delta_wb = sym('delta_wb', [3 1]);
    alpha_sym = sym('angular_acceleration', [3 1]);
    
    dt_sym = sym('dt', 1);
    acceleration = sym('unbiased_accelerometer', [3 1]);

    DCM = quatfromdcm(quat_sym);

    delta_q_sym = symbolic_quatexp(delta_theta.' / 2);
    %delta_q_sym = [1, delta_theta.' / 2]; % first order approximation

    quat_sym_t = symbolic_quatmultiply(delta_q_sym, quat_sym);
    DCM_sym_t = quatfromdcm(quat_sym_t);
    gyro_bias_sym_t = gyro_bias_sym + delta_wb;
    accelerometer_bias_sym_t = accelerometer_bias_sym + delta_ab;
    Xe_sym_t = Xe_sym + delta_Xe;
    Ve_sym_t = Ve_sym + delta_Ve;
    
    a_meas = accelerometer_true_measurements_sym + accelerometer_bias_sym_t;
    a_meas_ground = a_meas + DCM_sym_t.' * [0;0;g];
    w_meas = angular_velocity_sym + gyro_bias_sym_t + DCM_sym_t.' * [0; 0; 7.29e-5];
    %max_range = 34.9066; % 2000dps in rps
    %w_meas = min(max_range, max(-max_range, w_meas)); %saturation
    barometer = sea_level_pressure_sym * ( 1 + 0.0065 * Xe_sym_t(3) / temperature ) ^ (-9.80665 * 0.0289644 / (8.3144598 * 0.0065));

    ground_measurement = [a_meas_ground; w_meas; barometer];
    measurement = [a_meas; w_meas; barometer];
    
    kalman_state_vector = [delta_Xe; delta_Ve; delta_theta; delta_ab; delta_wb; angular_velocity_sym];

%     Hx_symbolic_on_ground = jacobian(ground_measurement, [Xe_sym; Ve_sym; quat_sym.'; accelerometer_bias_sym; gyro_bias_sym; angular_velocity_sym]);
    H_deltax_symbolic_on_ground = jacobian(ground_measurement, kalman_state_vector);
    H_deltax_symbolic = jacobian(measurement, kalman_state_vector);
%     Hx_symbolic = jacobian(measurement, [Xe_sym; Ve_sym; quat_sym.'; accelerometer_bias_sym; gyro_bias_sym; angular_velocity_sym]);
    % Hx_symbolic = jacobian(h(angular_velocity_sym, Xe_sym(3), gyro_bias_sym, accelerometer_bias_sym, accelerometer_true_measurements_sym, q2dcm_sym(quat(1), quat(2), quat(3), quat(4))), [Xe_sym; Ve_sym; quat; accelerometer_bias_sym; gyro_bias_sym; angular_velocity_sym]);
    
    global Hx_fun_ground Hx_fun
    
    Hx_fun_ground = matlabFunction(H_deltax_symbolic_on_ground);
    Hx_fun = matlabFunction(H_deltax_symbolic);
    
    fx = [delta_Xe + delta_Ve * dt_sym; delta_Ve + dt_sym * (-skew(DCM * acceleration) * delta_theta - DCM * delta_ab ); delta_theta - DCM * dt_sym * delta_wb; delta_ab; delta_wb; angular_velocity_sym + alpha_sym * dt_sym];
    Fx_sym = jacobian(fx, kalman_state_vector);
    Fx_fun = matlabFunction(Fx_sym);
    
    is_barometer_being_used = 1;

    known_moment = [0;0;0];
    Mb = [0;0;0];
    known_body_Acceleration = [0;0;15];
    Body_Acceleration = [0;0;15];

    for t = tstart:dt:time
        if t >= 0
            known_moment = [0.1 * sin(t * 2*pi * 1); 0.1 * sin(t * 2*pi * 1); 0.0005 * sin(t * 2*pi/2)];
            Mb = known_moment + [0.05 * cos(t * 2*pi * 1/10); 0.05 * sin(t * 2*pi * 1/40); 0.00005 * sin(t * 2*pi * 1/20)]; %normrnd(0, [0.01; 0.01; 0.00002] .^ 0.5, 3, 1); % + [t/100; -t/200; normrnd(0, 5e-4)];

            known_body_Acceleration = [0; 0; max(20, 0)] + [Mb(1:2) / 0.4; 0];
            Body_Acceleration = known_body_Acceleration + normrnd(0, [1; 1; 4] .^ 0.5, 3, 1);
        end

        known_alpha = I \ (known_moment - cross(nominal_state(1:3), I*nominal_state(1:3)));

        alpha = I \ (Mb - cross(w,I*w));

        q = quatnormalize(quatmultiply(q, quatexp([0; dt/2 * (w + w + alpha * dt)/2].')));
        w = w + alpha * dt;

        DCM = quat2dcm(quatinv(q));

        GlobalAcceleration = DCM * Body_Acceleration + [0;0;-g];

        Ve = Ve + dt * GlobalAcceleration;
        Xe = Xe + dt * Ve + dt^2/2 * GlobalAcceleration;

        if Xe(3) <= 0
            if Ve(3) < 0
                Ve(3) = 0;
            end
            Xe(3) = 0;
        end

        ISM330DHCX_noise = deg2rad(0.06^2);
        measurement_covariance = [.1 .1 .1 ISM330DHCX_noise ISM330DHCX_noise ISM330DHCX_noise 1].';
        measurement_covariance = measurement_covariance(1:6 + is_barometer_being_used);
        y = h(w, Xe(3), [0; 0; 0], [0; 0; 0], Body_Acceleration, DCM, is_barometer_being_used, g, 99364, 273.16 + 15) + normrnd(0, measurement_covariance .^ 0.5, length(measurement_covariance), 1);

        [nominal_state, P] = update(y, known_alpha, known_body_Acceleration, nominal_state, P, dt, diag(measurement_covariance), is_barometer_being_used, g, Fx_fun);
        error_w(:, round((t - tstart)/dt + 1)) = quat2axang(quatmultiply(q, quatinv(nominal_state(4:7).')));
        store_nominal_state(:, round((t - tstart)/dt + 1)) = nominal_state;
        
        if nominal_state(13) <= 0 && t > 0
            time = t;
            break
        end
    end

    rad2deg(error_w(4, end))

    figure(1)
    plot(tstart:dt:time, rad2deg(error_w(4, :)))
%     figure(2)
%     plot(tstart:dt:time, store_nominal_state(11:13, :))
%     legend('x','y', 'z')
end

function Hx = measurement_jacobian(height, delta_height, delta_theta, qw, qx, qy, qz, sea_level_pressure, temperature, is_barometer_being_used)
    global Hx_fun_ground Hx_fun
    if height <= 0
        Hx = Hx_fun_ground(height, delta_height, delta_theta(1), delta_theta(2), delta_theta(3), qw, qx, qy, qz, sea_level_pressure, temperature);
    else
        Hx = Hx_fun(height, delta_height, delta_theta(1), delta_theta(2), delta_theta(3), qw, qx, qy, qz, sea_level_pressure, temperature);
    end
    Hx = Hx(1:6 + is_barometer_being_used, :); % remove last row if barometer not being used
end

function [nominal_state, P] = update(y, alpha, Acceleration, nominal_state, P, dt, measurement_covariance, is_barometer_being_used, g, Fx_fun) 
    % http://www.iri.upc.edu/people/jsola/JoanSola/objectes/notes/kinematics.pdf

    w = nominal_state(1:3);
    q = nominal_state(4:7);
    % Ve = nominal_state(8:10);
    % Xe = nominal_state(11:13);
    % ab = nominal_state(14:16);
    % wb = nominal_state(17:19);
    % Thrust = nominal_state(20)
    sea_level_pressure = 99364;

    %Nominal state 236
    nominal_state(1:3) =  w + alpha * dt;
    nominal_state(4:7) = quatnormalize(quatmultiply(q.', quatexp([0; dt/2 * (w + nominal_state(1:3))/2].') ));

    DCM = quat2dcm(quatinv(nominal_state(4:7).'));

    GlobalAcceleration = DCM * Acceleration + [0;0;-g];

    nominal_state(8:10) = nominal_state(8:10) + dt * GlobalAcceleration; 
    nominal_state(11:13) = nominal_state(11:13) + dt * nominal_state(8:10) + dt^2/2 * GlobalAcceleration;

    if nominal_state(13) <= 0 && GlobalAcceleration(3) < 0
        GlobalAcceleration(3) = 0;
        nominal_state(13) = 0;
        nominal_state(10) = 0;
    end

    % nominal_state(14:16) = ab;
    % nominal_state(17:19) = wb;

    qw = nominal_state(4);
    qx = nominal_state(5); 
    qy = nominal_state(6); 
    qz = nominal_state(7);
    
    error_state = [zeros(15, 1); nominal_state(1:3)];

    %310
    % Fx = [eye(3), eye(3) * dt, zeros(3, 9); ...
    %     zeros(3,3), eye(3), - skew(DCM * (GlobalAcceleration)) * dt, -DCM * dt, zeros(3,3); ...
    %     zeros(3,6), eye(3), zeros(3,3), -DCM * dt; ...
    %     zeros(6,9), eye(6)];

%     Fx = [eye(3), eye(3) * dt, zeros(3, 9), zeros(3); ...
%         zeros(3,3), eye(3),  - skew(DCM * (GlobalAcceleration)) * dt, -DCM * dt, zeros(3,3), zeros(3); ...
%         zeros(3,6), eye(3), zeros(3,3), -DCM * dt, zeros(3); ...
%         zeros(3,9), eye(3), zeros(3, 6); ...
%         zeros(3,12), eye(3), zeros(3); ...
%         zeros(3,15), eye(3)];
    Fx = Fx_fun(dt, qw, qx, qy, qz, Acceleration(1), Acceleration(2), Acceleration(3));

    % 311
    Fi = [zeros(3, 15); eye(15)];

    %260
    w_process = diag([0.1 0.1 0.01]);
    acceleration_process = diag([0.2 0.2 0.5]);
    ISM330DHCX_bias_drift = (deg2rad(3)/3600)^2;
    gyro_bias_process = diag([ISM330DHCX_bias_drift ISM330DHCX_bias_drift ISM330DHCX_bias_drift]);
    %gyro_bias_proces = zeros(3);
    accel_bias_process = diag([1e-10 1e-10 1e-10]);
    %accel_bias_process = zeros(3);

    %270
    Qi = dt * [dt * acceleration_process, zeros(3, 12); ...
        zeros(3), dt * w_process, zeros(3, 9); ...
        zeros(3,6), accel_bias_process, zeros(3, 6); ...
        zeros(3, 9), gyro_bias_process, zeros(3); ...
        zeros(3,12), w_process];

    %267
    % deltax_hat = Fx * deltax_hat

    %268
    P = Fx*P*Fx.' + Fi * Qi * Fi.';

    %epsilon = y_t - y_t|t-1, h(x_t|t-1) = w_t|t-1
    measurement = h(nominal_state(1:3), nominal_state(13), nominal_state(17:19), nominal_state(14:16), Acceleration, DCM, is_barometer_being_used, g, sea_level_pressure, 273.16 + 15);

    % gyro_bias = sym('gyro_bias', [3 1]);
    % accelerometer_bias = sym('accelerometer_bias', [3 1]);
    % accelerometer_true_measurements = sym('accelerometer_true_measurements', [3 1]);
    % angular_velocity = sym('angular_velocity', [3 1]);
    % Xe = sym('Xe', [3 1]);
    % Ve = sym('Ve', [3 1]);
    % q = sym('q', [4 1], 'real');
    % syms q2dcm(qin);
    % qin = sym('qin', [1 4], 'real');
    % 
    % q2dcm(qin) = [qin(:,1).^2 + qin(:,2).^2 - qin(:,3).^2 - qin(:,4).^2 , 2.*(qin(:,2).*qin(:,3) + qin(:,1).*qin(:,4)), 2.*(qin(:,2).*qin(:,4) - qin(:,1).*qin(:,3)) ; 2.*(qin(:,2).*qin(:,3) - qin(:,1).*qin(:,4)), qin(:,1).^2 - qin(:,2).^2 + qin(:,3).^2 - qin(:,4).^2, 2.*(qin(:,3).*qin(:,4) + qin(:,1).*qin(:,2)) ; 2.*(qin(:,2).*qin(:,4) + qin(:,1).*qin(:,3)), 2.*(qin(:,3).*qin(:,4) - qin(:,1).*qin(:,2)), qin(:,1).^2 - qin(:,2).^2 - qin(:,3).^2 + qin(:,4).^2];
    % Hx_symbolic = jacobian([accelerometer_true_measurements + accelerometer_bias; angular_velocity + gyro_bias + q2dcm(q(1), q(2), q(3), q(4))' * [0; 0; 7.29e-5]; Xe(3)], [Xe; Ve; q; accelerometer_bias; gyro_bias; angular_velocity]);
    % Hx = double(subs(Hx_symbolic, quat, nominal_state(4:7)));

    % Hx = [zeros(3,10), eye(3), zeros(3, 6); ...
    %     zeros(3,13), eye(3), eye(3); ...
    %     0, 0, 1, zeros(1, 16)];

    % Hx = [zeros(3,10), zeros(3), zeros(3, 6); ...
    %     zeros(3,13), zeros(3), eye(3); ...
    %     0, 0, 1, zeros(1, 16)];

%     %312c
%     Qdeltatheta = 1/2 * [-qx -qy -qz; qw qz -qy; -qz qw qx; qy -qx qw];
%     
%     %279
%     Xdeltax = [eye(6), zeros(6,12); ...
%         zeros(4,6), Qdeltatheta, zeros(4,9); ...
%         zeros(9), eye(9)];
    H = measurement_jacobian(nominal_state(13), error_state(3), error_state(7:9), qw, qx, qy, qz, sea_level_pressure, 273.16 + 15, is_barometer_being_used);

    %272
    V = measurement_covariance;

%     %277
%     H = Hx * Xdeltax;

    %273-275
    K = (P * H.') / (H * P * H.' + V);
    error_state = error_state + K * (y - measurement);

    %Footnote 26
    %P = (eye(15) - K * H) * P;
    P = (eye(length(P)) - K * H) * P * (eye(length(P)) - K * H).' + K * V * K.';

    %313
    nominal_state(1:3) = error_state(16:18);
    nominal_state(11:13) = nominal_state(11:13) + error_state(1:3);
    nominal_state(8:10) = nominal_state(8:10) + error_state(4:6); 
    nominal_state(4:7) = quatnormalize( quatmultiply( quatexp([0; error_state(7:9) / 2].') , nominal_state(4:7).' ) );
    nominal_state(14:16) = nominal_state(14:16) + error_state(10:12);
    nominal_state(17:19) = nominal_state(17:19) + error_state(13:15);

    %315
    G = [eye(6), zeros(6,12); ...
        zeros(3,6), eye(3) + skew(error_state(7:9)/2), zeros(3,9); ...
        zeros(9), eye(9)];
    P = G*P*G.';
end

function y = skew(x)
    y = [0 -x(3) x(2) ; x(3) 0 -x(1) ; -x(2) x(1) 0 ];
end

function dcm = quatfromdcm(qin)
% assumes unit quaternion
    dcm = [qin(:,1).^2 + qin(:,2).^2 - qin(:,3).^2 - qin(:,4).^2 , 2.*(qin(:,2).*qin(:,3) + qin(:,1).*qin(:,4)), 2.*(qin(:,2).*qin(:,4) - qin(:,1).*qin(:,3)) ; ...
    2.*(qin(:,2).*qin(:,3) - qin(:,1).*qin(:,4)), qin(:,1).^2 - qin(:,2).^2 + qin(:,3).^2 - qin(:,4).^2, 2.*(qin(:,3).*qin(:,4) + qin(:,1).*qin(:,2)) ; ...
    2.*(qin(:,2).*qin(:,4) + qin(:,1).*qin(:,3)), 2.*(qin(:,3).*qin(:,4) - qin(:,1).*qin(:,2)), qin(:,1).^2 - qin(:,2).^2 - qin(:,3).^2 + qin(:,4).^2].';
end

function y = h(angular_velocity, height, gyro_bias, accelerometer_bias, accelerometer_true_measurements, DCM, is_barometer_being_used, g, sea_level_pressure, temperature)
    a_meas = accelerometer_true_measurements + accelerometer_bias;
    if height <= 0
        a_meas = a_meas + DCM.' * [0;0;g];
    end
    w_meas = angular_velocity + gyro_bias + DCM.' * [0; 0; 7.29e-5];
    max_range = 34.9066; % 2000dps in rps
    w_meas = min(max_range, max(-max_range, w_meas)); %saturation
    %w_meas = angular_velocity + gyro_bias;
    %barometer = height;
    % height = ((99364/Pressure)^(8.3144598*0.0065 / (9.80665*0.0289644)) - 1) * (273.16+temperature)/0.0065
    %barometer_temperature = temperature - 273.16;
    barometer = sea_level_pressure * ( 1 + 0.0065 * height / temperature ) ^ (-9.80665 * 0.0289644 / (8.3144598 * 0.0065));
    % temperature = temp - 273.16
    if is_barometer_being_used
        y = [a_meas; w_meas; barometer];
    else 
        y = [a_meas; w_meas];
    end
end

function quat = symbolic_quatmultiply(q, r)
    vec = [q(:,1).*r(:,2) q(:,1).*r(:,3) q(:,1).*r(:,4)] + ...
                 [r(:,1).*q(:,2) r(:,1).*q(:,3) r(:,1).*q(:,4)]+...
                 [ q(:,3).*r(:,4)-q(:,4).*r(:,3) ...
                   q(:,4).*r(:,2)-q(:,2).*r(:,4) ...
                   q(:,2).*r(:,3)-q(:,3).*r(:,2)];

    % Calculate scalar portion of quaternion product
    % scalar = s1*s2 - dot(v1,v2)
    scalar = q(:,1).*r(:,1) - q(:,2).*r(:,2) - ...
                 q(:,3).*r(:,3) - q(:,4).*r(:,4);
    quat = [scalar vec];
end

function quat = symbolic_quatexp(q)
    theta = sum(q .* q)^0.5;
    phi = sym('phi', 1);
    %44
    quat = subs(taylor([cos(phi), sinc(phi/pi) * q], 'Order', 6), phi, theta); %to prevent singularity when norm is 0 taylor series of sinx/x is used, taylor series of cos is also used as jacobian ends up dividing by q which leads to singularity when q = [0 0 0]
end