t = 0:0.001:5;
initial = [0 0 0 0 2 1]; % inital position and inital velocity

[t,y] = ode45(@test, t, initial)

h = animatedline;
%axis([min(y(:,1)), max(y(:,1)), min(y(:,2)), max(y(:,2)), min(y(:,3)), max(y(:,3))]);
view(3);

for k = 1:length(t)
    addpoints(h,y(k,1), y(k,2), y(k,3));
    drawnow
end

function y = test(t,x)
    y = zeros(6,1);
    y(1:3) = x(4:6);
    y(4) = y(1);
    y(5) = y(2) / 10;
    y(6) = y(3) / 5;
end