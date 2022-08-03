using TVC, DifferentialEquations, Plots, LinearAlgebra

using Parameters

@with_kw struct ODEParameters{R, V}
    veh::RocketParameters = RocketParameters()
    atmos::Atmosphere = Atmosphere(g = height -> [0; 0; 0])
    traj::RocketTrajectoryParameters = RocketTrajectoryParameters()

    Aero::Bool = false
    wind::V = zeros(3)

    ground::Bool = false

    MotorIgnitionTime::R = 0.0
    Control = (force=[1; 3; 5], torque=[0.25; 0.5; 0.75])
end

#   CASE 1
#   ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡

veh = RocketParameters(BurnTime=10., InertiaTensor=I(3), Mass=t -> 10)

x₀ = [[0; 0; 0]; [0; 0; 0]; [1; 0; 0; 0]; [0; 0; 0]][[veh.id_r; veh.id_v; veh.id_quat; veh.id_ω]];

p = ODEParameters(veh=veh)

#   Solve ODE
#   ≡≡≡≡≡≡≡≡≡≡≡

tspan = (0, veh.BurnTime)

prob = ODEProblem(f!, x₀, tspan, p)

sol = DifferentialEquations.solve(prob, reltol=1e-8, abstol=1e-8)

println(sol.u[end][veh.id_r])
println(TVC.Utils.to_matrix(sol.u[end][veh.id_quat]))
println(sol.u[end][veh.id_ω])

#   CASE 2
#   ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡

veh = RocketParameters(BurnTime=15., InertiaTensor=[0.7 0 0.7; 0 1 0; 0.7 0 -0.7], Mass=t -> 0.26)

x₀ = [[0; 0; 0]; [0; 0; 0]; [1; 0; 0; 0]; [0; 0; 0]][[veh.id_r; veh.id_v; veh.id_quat; veh.id_ω]];

p = ODEParameters(veh=veh, Control=(force=zeros(3), torque=[0; 0; 0.25]))

#   Solve ODE
#   ≡≡≡≡≡≡≡≡≡≡≡

tspan = (0, veh.BurnTime)

prob = ODEProblem(f!, x₀, tspan, p)

condition(x,t,integrator) = x[3] # when height is zero halt integration
cb = ContinuousCallback(condition, nothing, terminate!) # when going upwards do nothing

sol = DifferentialEquations.solve(prob, reltol=1e-8, abstol=1e-8)

println(sol.u[end][veh.id_r])
println(TVC.Utils.to_matrix(sol.u[end][veh.id_quat]))
println(sol.u[end][veh.id_ω])

#   CASE 3
#   ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡

veh = RocketParameters(BurnTime=15., InertiaTensor=[0.7 0 0.7; 0 1 0; 0.7 0 -0.7], Mass=t -> 0.26)

x₀ = [[0; 0; 0]; [0; 0; 0]; [1; 0; 0; 0]; [0; 10; 0]][[veh.id_r; veh.id_v; veh.id_quat; veh.id_ω]];


p = ODEParameters(veh=veh, Aero=false, 
    Control=(x, p, t) -> begin
        if t ≤ 5
            return (force=zeros(3), torque=[0.25; 0; -0.1])
        else 
            return (force=zeros(3), torque=zeros(3))
        end
    end)

#   Solve ODE
#   ≡≡≡≡≡≡≡≡≡≡≡

tspan = (0, veh.BurnTime)

prob = ODEProblem(f!, x₀, tspan, p)

sol = DifferentialEquations.solve(prob, reltol=1e-8, abstol=1e-8)

println(sol.u[end][veh.id_r])
println(TVC.Utils.to_matrix(sol.u[end][veh.id_quat]))
println(sol.u[end][veh.id_ω])