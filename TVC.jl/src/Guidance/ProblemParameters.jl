using LinearAlgebra, Parameters, StaticArrays
using ..Utils

export RocketTrajectoryParameters, RocketProblem

# ..:: Data structures ::..

""" Trajectory parameters. """
@with_kw struct RocketTrajectoryParameters{R<:Real} @deftype Union{SVector{3, R}, Vector{R}}
    r0 = @SVector Float64[20; -4; 30] # Initial Position
    v0 = @SVector Float64[4; -3; 0] # Initial Velocity
    q0::Union{SVector{4, R}, Vector{R}} = @SVector Float64[1; 0; 0; 0] # Initial Quaternion
    ω0 = @SVector zeros(3) # Initial Angular Velocity
    T0 = @SVector Float64[0; 0; 1] # Initial Thrust Vector
    Ṫ0 = @SVector zeros(3) # Initial Derivative of Thrust Vector
    # MotorFired::Base.RefValue{Bool} = Ref(false) # Has the motor been ignited yet, using ref value should allow mutability
    t0::R = 0.0 # time since motor was ignited.
    MotorFired::Bool = false

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