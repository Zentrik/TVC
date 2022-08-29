using Parameters, LinearAlgebra, StaticArrays
using ..Utils

export RocketParameters

@with_kw struct RocketParameters{R<:Real, M<:AbstractMatrix, I<:Int} @deftype R
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
    InertiaTensor::M = Diagonal([0.0826975856, 0.0826975856, 2.4778e-04]) # in body principal axis basis.

    Mass::Any = Mass
    Thrust::Any = Thrust
    Acceleration::Any = Acceleration
    CG::Any = CG
    COTToCG::Any = t -> Total_Length - CG(t)
    MomentArm::Any = t -> [0; 0; -COTToCG(t)]

    id_r::UnitRange{I} = 1:3
    id_v::UnitRange{I} = 4:6
    id_quat::UnitRange{I} = 7:10
    id_ω::UnitRange{I} = 11:13
    id_T::UnitRange{I} = 14:16
    id_Ṫ::UnitRange{I} = 17:19

    id_T̈::UnitRange{I} = 1:3
    id_roll::UnitRange{I} = 4:4

    id_tcoast::I = 1

    # nx::I = 19
    # nu::I = 4
    # np::I = 1
end