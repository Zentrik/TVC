inital_conditions = [pi/2 0]; %90 degrees and 0 angular velocity.
t = 0:0.001:5;

[t,y] = ode45( @threetest, t, inital_conditions);
y(:,3) = -9.81.*sin(y(:,1))./1; %alpha

close all; figure(); hold on
plot(t,y(:,1), 'b'); %theta
plot(t, y(:,2), 'r'); %omega
plot(t, y(:,3), 'g'); %alpha
xlabel('t');

function y = threetest(t,x)
    omega = x(2);
    alpha = -9.81 * sin(x(1))/1 -2*omega;
    y = [omega; alpha];
end
