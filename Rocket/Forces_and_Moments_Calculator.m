function [Forces, Moments]  = forces(Mach, Airspeed, AoA, Density, Temperature, DynamicPressure, CenterOfMassX, Thrust, body_angular_velocity, Roll_Torque, Xe)
    persistent ThetaRotation
    if isempty(ThetaRotation)
         ThetaRotation = eye(3);
    end
    
    persistent Croll_p
    if isempty(Croll_p)
        Croll_p = 0.03 * (rand - 0.5);
    end
    %% Rocket Shape Constants
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
   
    
    COT_to_COM = .95 - CenterOfMassX; %get actual value and change over time as solid motor burns
    
    %% General Constants
    K = 1.1;
   
    %% Coefficients
    angularrate_local = ThetaRotation \ body_angular_velocity;
    [CN, Cpitch] = CalculatePitchAndNormalCDs(Mach, K, Reference_Diameter, AoA, Airspeed, angularrate_local(2), NoseCone_Length, BodyTube_Length, Reference_Area, CenterOfMassX, PitchCenterX, NoseCone_PlanformArea, BodyTube_PlanformArea, NoseCone_CP, BodyTube_CP);

    Cbase = .12 + .13*Mach^2;
    CFrictionDrag = CalculateFrictionDrag(Temperature, Density, Reference_Radius, Reference_Area, SurfaceRoughness, NoseCone_Length, BodyTube_Length, Airspeed, Mach, NoseCone_WetArea, BodyTube_WetArea);
    Cpressure = PressureDragForce(sinphi);
    
    Cdrag = Cbase + CFrictionDrag + Cpressure;
    Caxial = AxialDragMul(AoA) * Cdrag;

    Cyaw = 0;
    Cside = 0;
    
     %% Random Noise to pitch and yaw coefficients in the range of +- 5e-4
    Cpitch = Cpitch + .1 * (rand - 0.5);
    Cyaw = Cyaw + .1 * (rand - 0.5);
    Croll = Croll_p + .0001 * (rand - 0.5);
    %% Forces and Moments
    lengtha = norm(Airspeed(1:2));
    
    cosa = Airspeed(1) / lengtha;
    sina = Airspeed(2) / lengtha;

    ThetaRotation = [cosa, -sina, 0; sina, cosa, 0; 0, 0, 1];
    
    Thrust_moment = COT_to_COM * [Thrust(1:2); 0] + [0;0; Roll_Torque];
    
    Forces = ThetaRotation * DynamicPressure * Reference_Area * -[CN; Cside; Caxial] + Thrust;
    Moments = ThetaRotation * (DynamicPressure * Reference_Area * Reference_Diameter * [-Cyaw; Cpitch; Croll]) + Thrust_moment;
    
    if Xe(3) <= 0
        Forces = [0;0;Thrust(3)];
        Moments = [0; 0; 0];
    end
    
end
