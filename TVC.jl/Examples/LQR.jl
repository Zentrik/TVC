using TVC, DifferentialEquations, Plots, LinearAlgebra, Parameters, ForwardDiff

#   Specify Parameters
#   ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡

veh = RocketParameters()
atmos = Atmosphere()

x₀ = [[0; 0; 0]; [0; 0; 0]; [1; 0; 0; 0]; [0; 0; 0]][[veh.id_r; veh.id_v; veh.id_quat; veh.id_ω]];

@with_kw struct ODEParameters{R, V}
    veh::RocketParameters = RocketParameters()
    atmos::Atmosphere = Atmosphere()

    Aero::Bool = false
    wind::V = zeros(3)

    ground::Bool = true

    MotorIgnitionTime::R = 0.0 # How long till motor ignites
    Control = (x, p, t) -> (force=[0; 0; veh.Thrust(t)], torque=zeros(3))
end

#   Continuous Actuator Model
#   ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡

function Actuator(x, p, t, desired_torque) # Model of TVC
    veh = p.veh
    tₘ = motorTime(t, p.MotorIgnitionTime)

    Thrustxy = [0 1; -1 0] * desired_torque[1:2] / veh.MomentArm(tₘ)[3];
    # r \times Thrust = Mb, if r = [0; 0; z] then Thrust = (Mb(2) / z, - Mb(1) / z, 0), then we add in z component of Thrust
    # (Mb(2) / z, - Mb(1) / z) = [0 1; -1 0] * Mb / z
    # See StateSpaceController.m for context

    if veh.Thrust(tₘ) ≈ 0
        Thrust = zeros(3)
    else
        if norm(Thrustxy) < veh.Thrust(tₘ) * sind(5)
            Thrustz = sqrt(veh.Thrust(tₘ)^2 - norm(Thrustxy)^2);
            Thrust = [Thrustxy; Thrustz]
        else 
            Thrustxy = normalize(Thrustxy) * veh.Thrust(tₘ) * sind(5)
            Thrustz = veh.Thrust(tₘ) * cosd(5)
            Thrust = [Thrustxy; Thrustz]
        end
    end

    roll = clamp(desired_torque[3], -.1, .1)
    torque = veh.MomentArm(tₘ) × Thrust + [0; 0; roll]
    
    return (force=Thrust, torque=torque)
end

#   LQR Controller
#   ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡

#   Specify Linearisation Points
#   ––––––––––––––––––––––––––––––

x̄ = x₀
ū = zeros(3)
t̄ = 1. # gives roughly mean thrust.

LQRindices = [veh.id_r[1:2]; veh.id_v[1:2]; veh.id_quat[2:4]; veh.id_ω] # states for LQR

p̄(u) = (veh=veh, atmos=atmos, traj=RocketTrajectoryParameters(), Aero=false, wind=zeros(3), ground=false, MotorIgnitionTime=0.0, Control=(x, p, t) -> Actuator(x, p, t, u)); # ground false to avoid the discontinuity, if touching ground, the forces and torques will be 0.

A = ForwardDiff.jacobian((dx, x) -> f!(dx, x, p̄(ū), t̄), zeros(13), x̄)
B = ForwardDiff.jacobian((dx, u) -> f!(dx, x̄, p̄(u), t̄), zeros(13), ū)

A = A[LQRindices, LQRindices] # janky way to ignore z component of r, v and real component of quat.
B = B[LQRindices, :]

using ControlSystems

Q = Diagonal([1,1,1, 2,2,2, 1,1, 1,1])
R = Diagonal([0.3, 0.3, 0.47])
K = lqr(Continuous, A, B, Q, R)

function control(x, p, t)
    return - K * x[LQRindices]
end

#   Solve ODE
#   ≡≡≡≡≡≡≡≡≡≡≡

tspan = (p.MotorIgnitionTime, p.MotorIgnitionTime + veh.BurnTime)

wind = randn(3) * 2
p = ODEParameters(veh=veh, atmos=atmos, Aero=true, wind=wind, Control=(x, p, t) -> Actuator(x, p, t, control(x, p, t)))

prob = ODEProblem(f!, x₀, tspan, p)

condition(x,t,integrator) = x[3] # when height is zero halt integration
cb = ContinuousCallback(condition, nothing, terminate!) # when going upwards do nothing

sol = DifferentialEquations.solve(prob, callback=cb)

# plot(sol, vars=veh.id_r)
# plot(sol, vars=veh.id_quat)
# plot(sol, vars=veh.id_ω)

deviationAngle = map(q -> acosd(TVC.Utils.rotate(q, [0; 0; 1]) ⋅ [0; 0; 1]), eachcol(sol[veh.id_quat, :])) # In degrees

plot(sol.t, deviationAngle, label="LQR")

# p = ODEParameters(veh=veh, atmos=atmos, Aero=true, wind=wind)

# prob = ODEProblem(f!, x₀, tspan, p)

# sol = DifferentialEquations.solve(prob, reltol=1e-8, abstol=1e-8, callback=cb)

# plot!(sol.t, map(q -> acosd(TVC.Utils.rotate(q, [0; 0; 1]) ⋅ [0; 0; 1]), eachcol(sol[veh.id_quat, :])), label="No LQR")