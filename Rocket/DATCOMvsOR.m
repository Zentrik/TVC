data = datcomimport('C:\Users\rag\Desktop\datcom.out');
data = data{1};

Reference_Diameter = 7.62e-2;
Reference_Radius = Reference_Diameter / 2;
Reference_Area = Reference_Radius^2 * pi;

NoseCone_Length = 15e-2;
BodyTube_Length = 80e-2; 
NoseCone_WetArea = pi * Reference_Radius * NoseCone_Length;          % pi*r*l 
BodyTube_WetArea = pi * Reference_Radius * 2 * BodyTube_Length;      % pi*d*l
NoseCone_PlanformArea = Reference_Diameter * NoseCone_Length / 2;
BodyTube_PlanformArea = Reference_Diameter * BodyTube_Length;

SurfaceRoughness = 20e-6;

sinphi = Reference_Radius / hypot(Reference_Radius, NoseCone_Length);
BodyTube_CP = 0.55;
NoseCone_CP = 0.1;
PitchCenterX = 0;
K = 1.1;

Mach = data.mach;
Airspeed = Mach * 340.064;
pitchrate = 0;
CenterOfMassX = 0.394;
AoA = (-pi/2:deg2rad(0.1):pi/2)';

NormalForceCD = zeros(1, length(AoA));
PitchMomentCD = zeros(1, length(AoA));

for i = 1:length(AoA)
    [NormalForceCD(i), PitchMomentCD(i)] = CalculatePitchAndNormalCDs(Mach, K, Reference_Diameter, AoA(i), Airspeed, pitchrate, NoseCone_Length, BodyTube_Length, Reference_Area, CenterOfMassX, PitchCenterX, NoseCone_PlanformArea, BodyTube_PlanformArea, NoseCone_CP, BodyTube_CP);  
end

[minValue,closestIndex] = min(abs(AoA - deg2rad(data.alpha)),[],1);

close all
subplot(2, 1, 1)
plot(data.alpha, data.cm(:, 1), rad2deg(AoA(closestIndex(1):closestIndex(end))), PitchMomentCD(closestIndex(1):closestIndex(end)))
legend('DATCOM', 'OpenRocket')
ylabel('CM')
xlabel('AoA')

subplot(2, 1, 2)
plot(data.alpha, data.cn(:, 1), rad2deg(AoA(closestIndex(1):closestIndex(end))), NormalForceCD(closestIndex(1):closestIndex(end)))
legend('DATCOM', 'OpenRocket')
ylabel('CN')
xlabel('AoA')