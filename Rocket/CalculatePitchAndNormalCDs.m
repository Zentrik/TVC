function [NormalForceCD, PitchMomentCD] = CalculatePitchAndNormalCDs(Mach, K, Reference_Diameter, AoA, Airspeed, pitchrate, NoseCone_Length, BodyTube_Length, Reference_Area, CenterOfMassX, PitchCenterX, NoseCone_PlanformArea, BodyTube_PlanformArea, NoseCone_CP, BodyTube_CP)
    cnmul = 1;
    if Mach < 0.05 && AoA > pi/4 
        cnmul = (Mach/0.05)^2;
    end
    NoseCone_CN = sin(AoA) * ( sin(AoA) * cnmul * K * NoseCone_PlanformArea / Reference_Area + 2 );
    BodyTube_CN = sin(AoA) * ( sin(AoA) * cnmul * K * BodyTube_PlanformArea / Reference_Area );
    NormalForceCD = BodyTube_CN + NoseCone_CN;

    CNCPxL = ( NoseCone_CN * NoseCone_CP + BodyTube_CN * BodyTube_CP) / Reference_Diameter;

    cacheDiameter = (NoseCone_PlanformArea + BodyTube_PlanformArea) / (NoseCone_Length + BodyTube_Length);

    mul = 3 * (PitchCenterX^4 + (NoseCone_Length + BodyTube_Length - PitchCenterX)^4 ) * (.275 * cacheDiameter / (Reference_Area * Reference_Diameter));
    PitchDampingMomentCD = min( ( mul * ( pitchrate /norm(Airspeed) )^2 ) , CNCPxL )  * sign(pitchrate);

    CNCmxL = NormalForceCD * CenterOfMassX / Reference_Diameter;
    PitchMomentCD = CNCPxL - CNCmxL - PitchDampingMomentCD;