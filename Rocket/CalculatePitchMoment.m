function PitchMomentCD = CalculatePitchMoment(K, Reference_Diameter, AoA, Airspeed, Pitchrate, CN, NoseCone_Length, BodyTube_Length, Reference_Area, PitchCenterX, NoseCone_PlanformArea, BodyTube_PlanformArea)
    NoseCone_CN = sin(AoA) * ( sin(AoA) * K * NoseCone_PlanformArea / Reference_Area + 2 );
    BodyTube_CN = sin(AoA) * ( sin(AoA) * K * BodyTube_PlanformArea / Reference_Area );
    CNCPxL = ( NoseCone_CN * NoseCone_CP + BodyTube_CN * BodyTube_CP) / Reference_Diameter;
    
    cacheDiameter = (NoseCone_PlanformArea + BodyTube_PlanformArea) / (NoseCone_Length + BodyTube_Length);
    
    mul = 3 * (PitchCenterX^4 + (NoseCone_Length + BodyTube_Length - PitchCenterX)^4 ) * (.275 * cacheDiameter / (Reference_Area * Reference_Diameter));
    PitchDampingMomentCD = min( ( mul * ( Pitchrate/norm(Airspeed) )^2 ) , CNCPxL )   * sign(Pitchrate);
    
    CNCmxL = CN * CenterOfMassx / Reference_Diameter;
    PitchMomentCD = CNCPxL - CNCmxL - PitchDampingMomentCD;