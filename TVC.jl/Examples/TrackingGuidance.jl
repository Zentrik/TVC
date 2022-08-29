using TVC, DifferentialEquations, SCPToolbox, Plots, Revise, LinearAlgebra

#   Specify Parameters
#   ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡

veh = RocketParameters()
atmos = Atmosphere()
traj = RocketTrajectoryParameters(r0=[20.; -4.; 30 + atmos.g(0)[3]*1.32487334729801753^2/2], v0=[4.; -3.; atmos.g(0)[3]*1.32487334729801753], t0=0.5, MotorFired=true);
# traj = RocketTrajectoryParameters();

mdl = RocketProblem(veh, atmos, traj)

servoΔt = 0.02 # Servo step rate
x0 = [traj.r0; traj.v0; traj.q0; traj.ω0][[veh.id_r; veh.id_v; veh.id_quat; veh.id_ω]];

#   Solve Optimal Trajectory
#   ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡

solution = solveProblem(mdl);
printSolution(solution)

#   Specify Controller
#   ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡

function control(x, p, t, servoΔt)
    veh = p.veh
    tₘ = motorTime(t, p.MotorIgnitionTime)
    sol = p.solution

    td = tₘ - tₘ % servoΔt

    if p.traj.t0 ≤ tₘ ≤ veh.BurnTime
        time = (td - p.traj.t0) / (veh.BurnTime - p.traj.t0)
        desired_tvc = sample(sol.xc, time)[veh.id_T]
        desired_roll = sample(sol.uc, time)[veh.id_roll]
    else
        desired_tvc = zeros(3)
        desired_roll = 0.
    end

    return Actuator(x, p, t, desired_tvc, desired_roll)
end

function control(x, p, t)
    veh = p.veh
    tₘ = motorTime(t, p.MotorIgnitionTime)
    sol = p.solution

    if p.traj.t0 ≤ tₘ ≤ veh.BurnTime
        time = (tₘ - p.traj.t0) / (veh.BurnTime - p.traj.t0)
        desired_tvc = sample(sol.xc, time)[veh.id_T]
        desired_roll = sample(sol.uc, time)[veh.id_roll]
    else
        desired_tvc = zeros(3)
        desired_roll = 0.
    end

    return Actuator(x, p, t, desired_tvc, desired_roll)
end

function Actuator(x, p, t, desired_tvc, desired_roll) # Model of TVC
    veh = p.veh
    tₘ = motorTime(t, p.MotorIgnitionTime)

    Thrust = desired_tvc * veh.Thrust(tₘ)
    torque = veh.MomentArm(tₘ) × Thrust + [0; 0; desired_roll]
    
    return (force=Thrust, torque=torque)
end

#   Simulate controller using optimal trajectory
#   ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡

using Parameters

# @with_kw struct ODEParameters{R, V}
#     veh::RocketParameters = RocketParameters()
#     atmos::Atmosphere = Atmosphere()
#     traj::RocketTrajectoryParameters = RocketTrajectoryParameters()

#     Aero::Bool = false
#     wind::V = zeros(3)

#     solution::SCPSolution

#     MotorIgnitionTime::R = 0.0
#     Control = zeros(3)
# end

@with_kw struct ODEParameters{R, S, V}
    veh::RocketParameters = RocketParameters()
    atmos::Atmosphere = Atmosphere()
    traj::RocketTrajectoryParameters = RocketTrajectoryParameters()

    Aero::Bool = false
    wind::V = zeros(3)

    solution::S

    MotorIgnitionTime::R = 0.0
    Control = (x, p, t) -> (force=[0; 0; veh.Thrust(t)], torque=zeros(3))
end

MotorIgnitionTime = solution.p[veh.id_tcoast]

if traj.MotorFired == true
    MotorIgnitionTime = 0.
end

tspan = (traj.t0, traj.t0 + veh.BurnTime + 10.0)

condition(x,t,integrator) = x[3] # when zero halt integration
cb = ContinuousCallback(condition, nothing, terminate!) # when going upwards do nothing

# saved_values = SavedValues(Float64, Float64)
# scb = SavingCallback((x, t, integrator) -> [integrator.p.Control(x, integrator.p, t).force; integrator.p.Control(x, integrator.p, t).torque[3]], saved_values, save_everystep = true)

# cbs = CallbackSet(cb, scb)

#   Continuous Controller
#   =======================

p = ODEParameters(veh=veh, atmos=atmos, traj=traj, solution=solution, MotorIgnitionTime=MotorIgnitionTime, Control=(x, p, t) -> control(x, p, t))

#   Discrete Controller
#   =====================

# p = ODEParameters(veh=veh, atmos=atmos, traj=traj, solution=solution, MotorIgnitionTime=MotorIgnitionTime, Control=(x, p, t) -> control(x, p, t, servoΔt))

# discretecb = PeriodicCallback((integrator) -> nothing, servoΔt);
# cbs = CallbackSet(cb, discretecb)

prob = ODEProblem(f!, x0, tspan, p)

sol = DifferentialEquations.solve(prob, reltol=1e-6, abstol=1e-6, callback=cb);

# println(sol.prob.p.MotorIgnitionTime)
println(sol.t[end] - sol.prob.p.MotorIgnitionTime)

using PlotlyJS
plotlyjs()
Plots.plot(sol, vars=mdl.veh.id_r)

t = LinRange(traj.t0, mdl.veh.BurnTime, 1000) .+ sol.prob.p.MotorIgnitionTime
Plots.plot!(t, mapreduce(k -> sample(solution.xc, k)[mdl.veh.id_r], hcat, LinRange(0, 1, 1000))')
Plots.plot!(sol.prob.p.MotorIgnitionTime .+ sol.prob.p.traj.t0 .+ mdl.veh.BurnTime * solution.td, solution.xd[mdl.veh.id_r, :]')
vline!([sol.prob.p.MotorIgnitionTime; sol.prob.p.MotorIgnitionTime + mdl.veh.BurnTime; sol.t[end]], linestyle=:dash)

Plots.plot(sol, vars=mdl.veh.id_v)
Plots.plot!(t, mapreduce(k -> sample(solution.xc, k)[mdl.veh.id_v], hcat, LinRange(0, 1, 1000))')

Plots.plot(sol, vars=mdl.veh.id_ω)
Plots.plot!(t, mapreduce(k -> sample(solution.xc, k)[mdl.veh.id_ω], hcat, LinRange(0, 1, 1000))')

Plots.plot(sol, vars=mdl.veh.id_quat)
Plots.plot!(t, mapreduce(k -> sample(solution.xc, k)[mdl.veh.id_quat], hcat, LinRange(0, 1, 1000))')

Plots.plot(t, mapreduce(k -> control(zeros(13), p, k).force / veh.Thrust(k), hcat, t)')
Plots.plot!(t, mapreduce(k -> sample(solution.xc, k)[mdl.veh.id_T], hcat, LinRange(0, 1, 1000))')

println(norm(sol.u[end][4:6]))
println(norm(sol(mdl.veh.BurnTime)[4:6]))

# plot(saved_values.t, reduce(hcat, saved_values.saveval)[3, :])
# plot!(saved_values.t, mapreduce(t -> control(x, p, t, servoΔt), hcat, saved_values.t)[1, :])
# plot!(saved_values.t, mapreduce(t -> control(motorTime(t, MotorIgnitionTime), veh, solution), hcat, saved_values.t)[1, :])

# Plots.plot(-1:0.1:5, reduce(hcat, sample.(Ref(solution.xc), [x / veh.BurnTime for x in -1:0.1:5]))[veh.id_T[1:3], :]')
# Plots.plot(0:0.01:5, reduce(hcat, [control(zeros(13), p, t).force for t in 0:0.01:5])')