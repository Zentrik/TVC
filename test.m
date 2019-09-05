z = [cos(pi/2/2), 0, 0, sin(pi/2/2)];
y = [cos(pi/8), 0, sin(pi/8), 0];
x = [cos(pi/16), sin(pi/16), 0, 0];
r = quatmultiply(z, y);
ori = quatmultiply(r, x);

theta = acos(ori(1));
x = atan(tan(theta) * ori(4)/sin(theta));
roll = [cos(x), 0, 0,sin(x)];
error = quatmultiply(quatinv(ori), roll);

quatmultiply(roll, quatinv(error));% is roll * rotation = orientation
quatmultiply(ori, error); % should be equal to roll, what the rocket position should end up as

w = [21/180*pi 21/180*pi  0];
W = [0 w];
q = [1 0 0 0];

%%
w_average = [2.5 1 0];
angle = norm(w_average);
axis_angle_of_rotation = [angle w_average/angle];
quaternion_axis_angle_of_rotation = [cos(angle/2) sin(angle/2)*axis_angle_of_rotation(2:4)];

%q_orientation = [1 0 0 0];
q_orientation = quatnormalize([cos(50) sin(50)*[.7071 0.7071 0]]);

if abs(q_orientation(1)) ~= 1 %is it rotated 0 degrees
    s = (1-q_orientation(1)^2)^0.5;
    angle_axis_orientation = [2*acos(q_orientation(1)) q_orientation(2:4)/s];
    y_error = angle_axis_orientation(1) * angle_axis_orientation(2);
    x_error = angle_axis_orientation(1) * angle_axis_orientation(3);
else 
    y_error = 0;
    x_error = 0;
end


q_new_orientation = quatmultiply(q_orientation, quaternion_axis_angle_of_rotation);
s = (1-q_new_orientation(1)^2)^0.5;
angle_axis_new_orientation = [2*acos(q_new_orientation(1)) q_new_orientation(2:4)/s];
y_new_error= angle_axis_new_orientation(1) * angle_axis_new_orientation(2);
x_new_error= angle_axis_new_orientation(1) * angle_axis_new_orientation(3);

x_error_change = x_new_error - y_error;
y_error_change = y_new_error - y_error;

%%
wx = 0.1;
wy = 1;
q_orientation = [1 0 0 0];
dt = .02;


w = [wy wx 0];
angle = norm(w) * dt/2;
axis_angle_of_rotation = [angle w/norm(w)];
quaternion_axis_angle_of_rotation = [cos(angle) sin(angle)*axis_angle_of_rotation(2:4)];

q_orientation = quatnormalize(quatmultiply(q_orientation, quaternion_axis_angle_of_rotation));

if abs(q_orientation(1)) ~= 1 %is it rotated 0 degrees
    s = (1-q_orientation(1)^2)^0.5;
    angle_axis_orientation = [2*acos(q_orientation(1)) q_orientation(2:4)/s]
    y_error = angle_axis_orientation(1) * angle_axis_orientation(2)
    x_error = angle_axis_orientation(1) * angle_axis_orientation(3)
else 
    y_error = 0
    x_error = 0
end