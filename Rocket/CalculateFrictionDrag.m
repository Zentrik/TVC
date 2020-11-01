function FrictionCd = CalculateFrictionDrag(Temperature, Density, Reference_Radius, Reference_Area, Roughness, NoseCone_Length, BodyTube_Length, Airspeed, Mach, NoseCone_WetArea, BodyTube_WetArea)
    
    KinematicViscosity = (0.0000037291 + 0.000000049944 * Temperature)/Density;
    Reynolds = norm(Airspeed) * (NoseCone_Length + BodyTube_Length) / KinematicViscosity;
    
    Cf = (1.5 * log(Reynolds) - 5.6)^-2;
    c1 = 1 - 0.1 * Mach^2;
    Cfc1 = Cf * c1;
    
    RoughnessCorrection = c1;
    NoseCone_RoughnessLimited = .032 * (Roughness / NoseCone_Length)^0.2 * RoughnessCorrection;
    BodyTube_RoughnessLimited = .032 * (Roughness / BodyTube_Length)^0.2 * RoughnessCorrection;
    
    NoseCone_Cf = max(Cfc1, NoseCone_RoughnessLimited);
    BodyTube_Cf = max(Cfc1, BodyTube_RoughnessLimited);
    
    BodyFriction = BodyTube_Cf * BodyTube_WetArea + NoseCone_Cf * NoseCone_WetArea;
    
    fB = (NoseCone_Length + BodyTube_Length + 0.0001) / Reference_Radius;
    correction = 1 + 1 / (2 * fB);
    
    FrictionCd = BodyFriction * correction / Reference_Area;