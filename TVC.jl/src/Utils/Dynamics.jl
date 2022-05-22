using LinearAlgebra
using ..Utils
import TVC: motorTime
using SCPToolbox

export f!

function f!(dx, x, p, t)
    f!(dx, x, p, t, p.Control)
end

function f!(dx, x, p, t, control) # to allow for jacobian to be calculated wrt control
    veh = p.veh
    atmos = p.atmos

    tₘ = motorTime(t, p.MotorIgnitionTime)

    if control isa Function
        ControlTorque = control(x, p, t).torque
        ControlForce = control(x, p, t).force
    else
        ControlTorque = control.torque
        ControlForce = control.force
    end

    r = x[veh.id_r]
    height = r[3]
    v = x[veh.id_v]
    quat = x[veh.id_quat]
    ω = x[veh.id_ω]
        
    Force = atmos.g(height) * veh.Mass(tₘ) + rotate(quat, ControlForce)
    Torque = ControlTorque

    if p.Aero
        Aero = CalculateAero(x, p, t)

        Force += Aero.force
        Torque += Aero.torque
    end

    if height <= 0 && Force[3] <= 0
        Force = zeros(3)
        Torque = zeros(3)
    end

    RigidBodyDynamics!(dx, x, p, tₘ, Force, Torque)
end

function RigidBodyDynamics!(dx, x, p, tₘ, Force, Torque)
    veh = p.veh

    r = x[veh.id_r]
    height = r[3]
    v = x[veh.id_v]
    quat = x[veh.id_quat]
    ω = x[veh.id_ω]

    dx[veh.id_r] = v
    dx[veh.id_v] = Force / veh.Mass(tₘ)
    dx[veh.id_quat] = 1/2 * quatL(quat) * [0; ω]
    dx[veh.id_ω] = veh.InertiaTensor \ (Torque - cross(ω, veh.InertiaTensor * ω))
end

function motorTime(t, MotorIgnitionTime)
    return t - MotorIgnitionTime
end