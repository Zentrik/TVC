using Parameters, LinearAlgebra
include("../../Guidance/src/Rocket_Acceleration.jl")

@with_kw struct RocketStruct{R<:Real} @deftype R
    Reference_Diameter = 7.62e-2
    Reference_Radius = Reference_Diameter / 2
    Reference_Area = Reference_Radius^2 * pi
    NoseCone_Length = 15e-2
    BodyTube_Length = 80e-2
    Total_Length = NoseCone_Length + BodyTube_Length

    NoseCone_WetArea = pi * Reference_Radius * NoseCone_Length           
    BodyTube_WetArea =  pi * Reference_Radius * 2 * BodyTube_Length 
    NoseCone_PlanformArea = Reference_Diameter * NoseCone_Length / 2
    BodyTube_PlanformArea = Reference_Diameter * BodyTube_Length

    SurfaceRoughness = 20e-6

    sinphi = Reference_Radius / hypot(Reference_Radius, NoseCone_Length)
    BodyTube_CP = 0.55
    NoseCone_CP = 0.1
    PitchCenterX = 0.0
    K = 1.1

    BurnTime = 3.45
    InertiaTensor::AbstractMatrix{R} = Diagonal([0.0826975856, 0.0826975856, 2.4778e-04])

    Mass::Any = Mass
    Thrust::Any = Thrust
    Acceleration::Any = Acceleration
    CG::Any = CG
    COTToCG::Any = t -> Total_Length - CG(t)
end