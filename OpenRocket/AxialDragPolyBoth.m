sn = deg2rad(17);
p2 = deg2rad(90);
%% AxialPoly1
matrix = [0, 0, 0, 1; sn^3, sn^2, sn, 1; 0 ,0, 1, 0; 3*sn^2, 2*sn, 1, 0];
a = matrix \ [1;1.3;0;0];

%% AxialPoly2

matrix2 = [sn^4, sn^3, sn^2, sn, 1 ;p2^4, p2^3, p2^2, p2, 1 ;4*sn^3, 3*sn^2, 2*sn, 1, 0 ;4*p2^3, 3*p2^2, 2*p2, 1, 0 ;12*p2^2, 6*p2, 2, 0, 0];
a2 = matrix2 \ [1.3;0;0;0;0];

%%
Error = zeros(1, length(axialdragpolynomialmultiplier(:,1)));
MulList = zeros(1, length(axialdragpolynomialmultiplier(:,1)));
for i = 1:length(axialdragpolynomialmultiplier(:,1))
    AoA = deg2rad(axialdragpolynomialmultiplier(i,1));

    if AoA > p2
        AoA = pi - AoA;
    end

    if AoA < sn
        mul = dot(a, [AoA^3, AoA^2, AoA, 1]);
    else 
        mul = dot(a2, [AoA^4, AoA^3, AoA^2, AoA, 1]);
    end

    if deg2rad(axialdragpolynomialmultiplier(i,1)) >= p2
        mul = -mul;
    end
    
    if axialdragpolynomialmultiplier(i,2) ~= 0
        Error(i) = (axialdragpolynomialmultiplier(i,2) - mul)/axialdragpolynomialmultiplier(i,2);
    end
    
    MulList(i) = mul;
end

Error = 100 * Error;

close all;
figure(1);
plot(axialdragpolynomialmultiplier(:,1), Error);
xlabel('AoA')
ylabel('Error')

figure(2);
hold on
plot(axialdragpolynomialmultiplier(:,1), MulList, 'b');
plot(axialdragpolynomialmultiplier(:,1), axialdragpolynomialmultiplier(:,2), 'g');
hold off
xlabel('AoA')
legend('Matlab','OpenRocket')

mean(Error)
median(Error)
std(Error)