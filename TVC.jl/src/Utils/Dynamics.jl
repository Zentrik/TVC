using LinearAlgebra
using ..Utils
import TVC: motorTime
using SCPToolbox

export f!

function f!(dx, x, p, t)
    veh = p.veh
    atmos = p.atmos

    tₘ = motorTime(t, p.MotorIgnitionTime)

    if p.Control isa Function
        Thrust = p.Control(t)[1:3]
        roll = p.Control(t)[4]
    else
        Thrust = p.Control[1:3]
        roll = p.Control[4]
    end

    r = x[veh.id_r]
    height = r[3]
    v = x[veh.id_v]
    quat = x[veh.id_quat]
    ω = x[veh.id_ω]

    Acceleration = atmos.g(height) + rotate(quat, Thrust) / veh.Mass(tₘ)
    Torque = veh.MomentArm(tₘ) × Thrust + [0; 0; roll]

    RigidBodyDynamics!(dx, x, p, Acceleration, Torque)

    # CalculateAero!(dx[veh.id_v], dx[veh.id_ω], veh, p.atmos, p.wind, v, quat, ω, height, t)
end

function fTest!(dx, x, p, t)
    veh = p.veh
    atmos = p.atmos
    solution = p.solution

    tₘ = motorTime(t, p.MotorIgnitionTime)

    if p.Control isa Function
        T̂ = p.Control(t)
    else
        T̂ = p.Control[1:3]
    end

    r = x[veh.id_r]
    height = r[3]
    v = x[veh.id_v]
    quat = x[veh.id_quat]
    ω = x[veh.id_ω]

    Acceleration = atmos.g(height) + rotate(quat, T̂) * veh.Acceleration(tₘ)
    
    if 0. ≤ tₘ ≤ veh.BurnTime
        roll = sample(solution.uc, tₘ / veh.BurnTime)[veh.id_roll]
    else
        roll = 0.
    end

    Torque = veh.MomentArm(tₘ) × T̂ * veh.Thrust(tₘ) + [0; 0; roll]

    RigidBodyDynamics!(dx, x, p, Acceleration, Torque)

    # CalculateAero!(dx[veh.id_v], dx[veh.id_ω], veh, p.atmos, p.wind, v, quat, ω, height, t)
end

function RigidBodyDynamics!(dx, x, p, Acceleration, Torque)
    veh = p.veh
    atmos = p.atmos

    r = x[veh.id_r]
    height = r[3]
    v = x[veh.id_v]
    quat = x[veh.id_quat]
    ω = x[veh.id_ω]

    dx[veh.id_r] = v
    dx[veh.id_v] = Acceleration
    dx[veh.id_quat] = 1/2 * quatL(quat) * [0; ω]
    dx[veh.id_ω] = veh.InertiaTensor \ (Torque - cross(ω, veh.InertiaTensor * ω))
end

function motorTime(t, MotorIgnitionTime)
    return t - MotorIgnitionTime
end