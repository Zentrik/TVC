t = 0:0.001:5;
initial = [0 1]; % inital position and inital velocity

[t,y] = ode45(@test, t, initial);
y(:,3) = 5*y(:,2);  %calculates acceleration

close all; figure(); hold on
plot(t, y(:,1), 'b'); %position
plot(t, y(:,2), 'r'); %velocity
plot(t, y(:,3), 'g'); %acceleration
xlabel('t');
y(end,1)

function y = test(t,x)
    velocity = x(2);
    acceleration = velocity * 5;
    y = [velocity; acceleration]; %actually will be [position, acceleration], as ode integrates
end