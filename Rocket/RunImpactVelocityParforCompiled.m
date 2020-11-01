ThrustTable = xlsread('F15_Thrust','Sheet1');
mass = 1.0565;
ThrustAcceleration = griddedInterpolant([-1; ThrustTable(:,1); 100], [0; ThrustTable(:,2); 0] / mass);

s = 1:0.1:500;
v = 50:-0.1:-100;

data = ImpactVelocityParforCompile_mex(v, s, ThrustAcceleration(0:0.001:30), 0.001, 30, 3.45, 8);
% 
% % save ImpactVelocity data -v7.3
surf(v, s, data, 'EdgeColor', 'none');
xlabel('Initial Velocity (m/s)');
ylabel('Initial Height (m)');
zlabel('Maximum Impact Velocity (m/s)')
colormap(flipud(winter));

% p00 = 30.55; 
% p10 = 0.8126;
% p01 = 0.1106;
% p20 = 0.01797;
% p11 = 0.00196;
% p02 = 0.0007978;
% p30 = -0.0001494;
% p21 = -6.561e-05;
% p12 = -9.147e-06;
% p03 = -3.186e-06;
% p40 = -1.505e-06;
% p31 = -9.441e-08;
% p22 = -1.068e-09;
% p13 = -1.527e-09;
% p04 = 4.961e-09;
% p50 = 2.375e-09;
% p41 = 4.713e-09;
% p32 = 1.098e-09;
% p23 = 1.674e-10;
% p14 = 1.645e-11;
% p05 = -3.134e-12;
%           
% syms f(x,y)
% f(x,y) = p00 + p10*x + p01*y + p20*x^2 + p11*x*y + p02*y^2 + p30*x^3 + p21*x^2*y + p12*x*y^2 + p03*y^3 + p40*x^4 + p31*x^3*y + p22*x^2*y^2 + p13*x*y^3 + p04*y^4 + p50*x^5 + p41*x^4*y + p32*x^3*y^2 + p23*x^2*y^3 + p14*x*y^4 + p05*y^5;
% 
% numv = length(vz);
% nums = length(sz);
% numtotal = nums*numv;
% [vz0, sz0] = meshgrid(vz, sz);
% fz = c(vz0, sz0);

% surf(vz, sz, z, 'FaceColor', 'r', 'EdgeColor', 'none');
% hold on
% surf(vz, sz, fz, 'EdgeColor', 'none');
% xlabel('Initial Velocity (m/s)');
% ylabel('Initial Height (m)');
% fz = double(f(vz0, sz0));

% fz = zeros(nums, numv);
% 
% for j = 1:numtotal
%     fz(j) = double(f(sz(1 + mod(nums - 1 + j, nums)), vz(1 + floor((j - 1)/nums))));
% end

% function d = c(x, y)
%     p00 = 30.55; 
%     p10 = 0.8126;
%     p01 = 0.1106;
%     p20 = 0.01797;
%     p11 = 0.00196;
%     p02 = 0.0007978;
%     p30 = -0.0001494;
%     p21 = -6.561e-05;
%     p12 = -9.147e-06;
%     p03 = -3.186e-06;
%     p40 = -1.505e-06;
%     p31 = -9.441e-08;
%     p22 = -1.068e-09;
%     p13 = -1.527e-09;
%     p04 = 4.961e-09;
%     p50 = 2.375e-09;
%     p41 = 4.713e-09;
%     p32 = 1.098e-09;
%     p23 = 1.674e-10;
%     p14 = 1.645e-11;
%     p05 = -3.134e-12;
%     d = p00 + p10.*x + p01.*y + p20.*x.^2 + p11.*x.*y + p02.*y.^2 + p30.*x.^3 + p21.*x.^2.*y + p12.*x.*y.^2 + p03.*y.^3 + p40.*x.^4 + p31.*x.^3.*y + p22.*x.^2.*y.^2 + p13.*x.*y.^3 + p04.*y.^4 + p50.*x.^5 + p41.*x.^4.*y + p32.*x.^3.*y.^2 + p23.*x.^2.*y.^3 + p14.*x.*y.^4 + p05.*y.^5;
% end