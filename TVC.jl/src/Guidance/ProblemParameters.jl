using LinearAlgebra, Parameters, StaticArrays
using ..Utils
import TVC: motorTime

export RocketTrajectoryParameters, RocketProblem, horizon

# ..:: Data structures ::..

""" Trajectory parameters. """
@with_kw struct RocketTrajectoryParameters{R<:Real} @deftype Union{SVector{3, R}, Vector{R}}
    r0 = @SVector Float64[20; -4; 30] # Initial Position
    v0 = @SVector Float64[4; -3; 0] # Initial Velocity
    q0::Union{SVector{4, R}, Vector{R}} = @SVector Float64[1; 0; 0; 0] # Initial Quaternion
    ω0 = @SVector zeros(3) # Initial Angular Velocity
    T0 = Float64[0; 0; 1] #@SVector Float64[0; 0; 1] # Initial Thrust Vector
    Ṫ0 = zeros(3) #@SVector zeros(3) # Initial Derivative of Thrust Vector
    # MotorFired::Base.RefValue{Bool} = Ref(false) # Has the motor been ignited yet, using ref value should allow mutability
    t0::R = 0. # time since motor started firing (if MotorFired is true). # confusing that traj and odeparamaters both have t0
    MotorFired::Bool = false
    ExpectedIgnitionTime::R = 0. #(0.7 - 0.419), # varies from 0.605 - 0.419 to 0.8 - 0.419
    # Don't change from 0
    
    PreviousTrajectoryState::ContinuousTimeTrajectory = ContinuousTimeTrajectory([0], [0], :linear)
    PreviousTrajectoryInput::ContinuousTimeTrajectory = ContinuousTimeTrajectory([0], [0], :linear)
    PreviousTrajectoryCurrentTime::R = 0. # Current time relative to previous trajectory, i.e. sample(, current time) gives the previous trajectory's solution at the current time.
    PreviousTrajectoryP:: R = 0. # Best estimate for t_coast relative to current time.
    PreviousTrajectoryTLand::R = 3.45 # Best estimate for the touchdown motor time, i.e. the previous solution's p[veh.id_tland]. Defaults to the F15 burn time.
    UsePreviousTrajectory::Bool = false

    rN = @SVector zeros(3) # Final Position
    vN = @SVector zeros(3) # Final Velocity
    qN::Union{SVector{4, R}, Vector{R}} = @SVector Float64[1; 0; 0; 0] # Final Quaternion
    ωN = @SVector zeros(3) # Final Angular Velocity
    TN = @SVector Float64[0; 0; 1] # Final Thrust Vector
    ṪN = @SVector zeros(3) # Final Derivative of Thrust Vector
end

""" Rocket trajectory optimization problem parameters all in one. """
@with_kw struct RocketProblem
    veh::RocketParameters = RocketParameters() # The parameters of the rocket
    atmos::Atmosphere = Atmosphere() # The environment
    traj::RocketTrajectoryParameters = RocketTrajectoryParameters() # The trajectory
end

"""
    horizon(p, mdl)

How much powered flight the trajectory covers: the motor time at touchdown,
which is a decision variable, less the motor time we are planning from.
"""
horizon(p, mdl::RocketProblem) = p[mdl.veh.id_tland] - mdl.traj.t0

"""
    motorTime(τ, p, mdl)

Motor time at normalised trajectory time `τ ∈ [0, 1]`.
"""
motorTime(τ, p, mdl::RocketProblem) = mdl.traj.t0 + τ * horizon(p, mdl)