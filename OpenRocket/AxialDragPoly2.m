%% AxialPoly1
matrix = [0, 0, 0, 1; (17/180*pi)^3, (17/180*pi)^2, 17/180*pi, 1; 0 ,0, 1, 0; 3*(17/180*pi)^2, 2*(17/180*pi), 1, 0];
a = matrix \ [1;1.3;0;0];

AoA = [0.065517914790615];
y_actual = [1.03743191861063];
error = zeros(1, length(AoA));
for i = 1:length(AoA)
    x = [AoA(i)^3, AoA(i)^2, AoA(i), 1];
    y = dot(a, x);
    error(i) = 100 * (y - y_actual(i))/y_actual(i);
end


%% AxialPoly2
matrix2 = [(17/180*pi)^4, (17/180*pi)^3, (17/180*pi)^2, (17/180*pi), 1 ;(pi/2)^4, (pi/2)^3, (pi/2)^2, pi/2, 1 ;4*(17/180*pi)^2, 3*(17/180*pi)^2, 2*(17/180*pi), 1, 0 ;4*(pi/2)^3, 3*(pi/2)^2, pi, 1, 0 ;12*(pi/2)^2, 3*pi, 2, 0, 0];
a2 = matrix2 \ [1.3;0;0;0;0];

AoA2 = [1.29835043055858,0.914884140602907,0.392262749385726];
y_actual2 = [0.042686764571741, 0.435538681478357, 1.26038696537678];
error2 = zeros(1, length(AoA2));
for i = 1:length(AoA2)
    x2 = [AoA2(i)^4,  AoA2(i)^3, AoA2(i)^2, AoA2(i), 1];
    y2 = dot(a2, x2);
    error2(i) = 100 * (y2 - y_actual2(i))/y_actual2(i);
end

error, error2
