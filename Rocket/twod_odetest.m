t = 0:0.001:5;
initial = [0 0 2 1]; % inital position and inital velocity

[t,y] = ode45(@test, t, initial);
y(:,5) = y(:,3) * 2;;
y(:,6) = y(:,4) *2;

%{
close all;
plot3(t, y(:,1), y(:,2), 'b', t, y(:,3), y(:,4), 'g', t, y(:,5), y(:,6),'r');
xlabel('t');
ylabel('x');
zlabel('y');
y(end,1:2)
%}

h = animatedline;
axis([min(y(:,1)), max(y(:,1)), min(y(:,2)), max(y(:,2))])

for k = 1:length(t)
    addpoints(h,y(k,1), y(k,2));
    drawnow
end

function y = test(t,x)
    v1 = x(3);
    v2 = x(4);
    a1 = v1 * 2;
    a2 = v2;
    y = [v1; v2; a1; a2]; %actually will be [position, acceleration], as ode integrates
end