struct RocketStruct
    Reference_Diameter
    Reference_Radius
    Reference_Area
    NoseCone_Length
    BodyTube_Length

    NoseCone_WetArea       
    BodyTube_WetArea
    NoseCone_PlanformArea
    BodyTube_PlanformArea 

    SurfaceRoughness

    sinphi
    BodyTube_CP
    NoseCone_CP
    PitchCenterX
    K

    BurnTime
    InertiaTensor

    Mass
    Thrust
    Acceleration
    CG
end   

function RocketContsruct(Reference_Diameter, NoseCone_Length, BodyTube_Length, SurfaceRoughness, BodyTube_CP, NoseCone_CP, PitchCenterX, K, Burn_Time, InertiaTensor, Mass, Thrust, Acceleration, CG)
    Reference_Radius = Reference_Diameter / 2;
    Reference_Area = Reference_Radius^2 * pi;

    NoseCone_WetArea = pi * Reference_Radius * NoseCone_Length;        
    BodyTube_WetArea = pi * Reference_Radius * 2 * BodyTube_Length;
    NoseCone_PlanformArea = Reference_Diameter * NoseCone_Length / 2;
    BodyTube_PlanformArea = Reference_Diameter * BodyTube_Length;

    sinphi = Reference_Radius / hypot(Reference_Radius, NoseCone_Length)

    RocketStruct(Reference_Diameter, Reference_Radius, Reference_Area, NoseCone_Length, BodyTube_Length, NoseCone_WetArea, BodyTube_WetArea, NoseCone_PlanformArea, BodyTube_PlanformArea, SurfaceRoughness, sinphi, BodyTube_CP, NoseCone_CP, PitchCenterX, K, Burn_Time, InertiaTensor, Mass, Thrust, Acceleration, CG)
end