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

    # A solid motor cannot throttle. Set this to true to let the guidance
    # problem choose ‖T‖ ∈ [0, 1] anyway, which is what it used to do.
    Throttleable::Bool = false

    # Set this to true to force touchdown to happen exactly at burnout, which is
    # what the guidance problem used to require. With no throttle that leaves it
    # with almost nothing to steer the terminal altitude constraint with, so by
    # default the touchdown time is a decision variable and the rocket is allowed
    # to fall ballistically once the motor is spent.
    FixedLandingTime::Bool = false

    # Weight on the running cost, in units of the terminal cost (m²/s²): holding
    # one input channel at its limit for the whole flight costs this much, so
    # 0.01 is worth 0.1 m/s of touchdown speed. Without it nothing in the
    # problem prefers a small input at all.
    #
    # It buys a smoother gimbal command — measured on the nominal solve,
    # max‖T̈‖ goes 0.1745 (saturated) at w = 0.01, to 0.0940 at w = 1, to 0.0269
    # at w = 100, against touchdown speeds of 0.524, 0.571 and 0.626 m/s. Pick
    # a point on that curve to taste.
    #
    # It does *not* fix the roll rate. max|u₄| is 1e-4 to 1e-3 at every weight
    # above, so that residual is a numerical floor, not something the optimiser
    # is choosing — see docs/mpc-feasibility.md.
    InputCostWeight = 0.01

    # Loose bound on ‖ω‖, a safety net rather than a design constraint. The
    # value that used to be here (commented out) was π/2, which is below rates
    # that legitimately show up in recorded flight states, so it would have
    # fought the initial condition rather than shaping the trajectory.
    MaxAngularVelocity = 2 * pi

    # How long after burnout the guidance will plan an unpowered fall for.
    # Touchdown is never allowed before BurnTime: thrust to weight is ~1.35, so
    # a rocket that reaches the ground under thrust takes off again.
    MaxBallisticTime = 2.0

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
    id_tland::I = 2 # motor time at touchdown, see FixedLandingTime

    # nx::I = 19
    # nu::I = 4
    # np::I = 1
end