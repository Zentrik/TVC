quat = [cos(pi/4), sin(pi/4)^2, 0 ,sin(pi/4)^2];
desired = [1 0 0 0];
desired_roll = [cos(deg2rad(5)) sin(deg2rad(5)) * [cos(pi/4) cos(pi/4) 0]];
b_d_roll = quatmultiply(quatmultiply(quatinv(quat), desired_roll), quat);

k = [0 0 1];
vector = quatrotate(quatinv(quat), k);
vector_roll = quatrotate(quatinv(desired_roll), k);

rotation = quatmultiply(desired, quatinv(quat));

theta = atan( rotation(4) / rotation(1) ); %half the roll angle
roll = [cos(theta), 0, 0, sin(theta)]; %extract roll part of quaternion
no_roll = quatmultiply(quatinv(roll), rotation);

test = quatnormalize([1 + dot(vector, k) cross(vector, k)]); %rotation from vector to k
inertial_roll = quatmultiply(quatmultiply(quat, roll), quatinv(quat));  %roll about body axis in quat form using inertial axis

for theta = 1:0.001:180
       roll = [cos(deg2rad(theta)) 0 0 sin(deg2rad(theta))];
       quatie = quatmultiply(quatmultiply(no_roll, roll), quat);
       if abs(quatie(4)) < x
          quatie(4)
          theta
       end
end
 
%working
desired = desired_roll;
rotation_b = quatmultiply(quatinv(quat), desired);

theta = atan( rotation(4) / rotation(1) ); %half the roll angle
roll = [cos(theta), 0, 0, sin(theta)]; %extract roll part of quaternion
noroll_b = quatmultiply(rotation_b, quatinv(roll));

final_quat = quatmultiply(quat, noroll_b);
quatrotate(quatinv(final_quat), k)
quatrotate(quatinv(desired), k)
