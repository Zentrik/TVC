%{
Cs = 0
function alpha = coeficients(rho, v, A_ref, L_ref, AoA)

Cm = Cm - CN * centerofmass / Length
end

function [a, alpha] = acceleration(rho, v_wind, A_ref, CN ,Cs, mass, Thrust, Ca, q_A2B, I, Length, Cm, Cyaw, Croll, velocity)
v_airflow = norm(v_wind - velocity); 
z = rho * v_airflow^2 * A_ref/2;

a_local = [-CN * z, -Cs * z, Thrust - Ca * z] / mass;
a = quatrotate(q_A2B, a_local);
a(3) = a(3) - 9.81;
velocity = integrate acceleraiotneaiosneaoirntseionratoinrasitenrieostnersteiarstnaersnteiarnstoarnteioasno

moment = z * Length * [-Cyaw, Cm, Croll];
alpha = [moment(1)/I(1), moment(2)/I(1), moment(3)/I(2)];
end
%}
data = xlsread('/home/rahul/Documents/school/Current/EPQ/Rocket/F10_Thrust.xls','Sheet1');

TimeData = data(:,1);
ThrustData = data(:,2);

t = 0:0.001:20;
initial = [0 0 0 0 0 0]; %position, veloctiy
A_ref = 0.456;
q_A2B = [1 0 0 0];
pressure = 101325;
wind_velocity = [1; -1; 0];

Options = odeset('Events', @myEvent);
[t,y] = ode45(@(t,y) ODE(t, y, TimeData, ThrustData), t, initial, Options);
AoA = getAoA(y, wind_velocity);

close all;

figure(1);
plot(t, AoA);
%{
h = animatedline;
xlabel('x');
ylabel('y');
zlabel('z');
view(3);

for k = 1:length(t)
    addpoints(h,y(k,1), y(k,2), y(k,3));
    drawnow
end
%}

figure(2);
plot3(y(:,1), y(:,2), y(:,3));
xlabel('x');
ylabel('y');
zlabel('z');
y(end,1:3)


%{
close all;
plot(t, y(:,3));
%}

function y = ODE(t, x, TimeData, ThrustData)
    y = zeros(6,1);
    wind_velocity = [1; -1; 0];
    rho = 101325/287.053/293.15;
    A_ref= 0.456;
    Cs = 0;
    Ca = 0.2;
    mass = .88;
    q_A2B = [1 0 0 0];
    Length = 0.80;
    
    ThrustForce = interp1(TimeData, ThrustData, t);
    if isnan(ThrustForce)
        ThrustForce = 0;
    end
    
    y(1:3) = x(4:6); %inital velocity
    airflow_velocity = wind_velocity - x(4:6); 
    airflow_local = quatrotate(quatinv(q_A2B), airflow_velocity.');
    mag_airflow_local = norm(airflow_local);
    
    len = hypot(airflow_local(1), airflow_local(2));
    if mag_airflow_local > 0.01
        AoA = acos(airflow_local(3)/len);
    else 
        AoA = 0;
    end
    v_drag = mag_airflow_local^2;
    
    %CN = sin(AoA) * (1.1 * 212016280 * pi * sin(AoA)/ 43428647 + 2);
    CN = 0.5;
    
    z = rho * v_drag * A_ref/2;
    a_local = [-CN * z, -Cs * z, ThrustForce - Ca * z] / mass;
    y(4:6) = a_local.';
    %acceleration = quatrotate(q_A2B, a_local).'
    if x(3) < 1.e-6 && ThrustForce <= 9.81
        y(6) = 0;
    else
        y(6) = y(6) - 9.81;
    end
    
    W = [ 0 angular_velocity]
    q_A2B_dot = 
    moment = z*Length * [-Cyaw Cm Croll];
    alpha = moment .* [I(1) I(1) I(3)];
    
    
end

function [value, isterminal, direction] = myEvent(T, Y)
    value = Y(3);
    isterminal = 1;   % Stop the ode once landed
    direction = -1;
end

function AoA = getAoA(y, wind_velocity)
    q_A2B = [1 0 0 0];
    AoA = zeros(size(y,1), 1);
    for k = 1:size(y,1)
        airflow_velocity = wind_velocity - y(k,4:6); 
        airflow_local = quatrotate(quatinv(q_A2B), airflow_velocity);

        len = hypot(airflow_local(1), airflow_local(2));
        mag_airflow_local = norm(airflow_local);
        if mag_airflow_local > 0.01
            AoA(k) = acos(airflow_local(3)/len);
        end
    end
end