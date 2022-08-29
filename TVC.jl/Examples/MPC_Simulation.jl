using TVC, DifferentialEquations, SCPToolbox, Plots, Revise, LinearAlgebra, Parameters

#   Specify Parameters
#   ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡

veh = RocketParameters()
atmos = Atmosphere()
traj = RocketTrajectoryParameters();

mdl = RocketProblem(veh, atmos, traj)

servoΔt = 0.02 # Servo step rate
x_0 = [traj.r0; traj.v0; traj.q0; traj.ω0][[veh.id_r; veh.id_v; veh.id_quat; veh.id_ω]];

#   Solve Optimal Trajectory
#   ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡

solution = solveProblem(mdl);

#   Specify Controller
#   ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡

function control(x, p, t, servoΔt)
    veh = p.veh
    tₘ = motorTime(t, p.MotorIgnitionTime)
    sol = p.solution[]

    td = tₘ - tₘ % servoΔt

    if 0 ≤ tₘ ≤ veh.BurnTime
        desired_tvc = sample(sol.xc, td / veh.BurnTime)[veh.id_T]
        desired_roll = sample(sol.uc, td / veh.BurnTime)[veh.id_roll]
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

#   MPC Controller
#   ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡

function mpc!(integrator)
    veh = integrator.p.veh
    atmos = integrator.p.atmos

    x = integrator.u

    r = x[veh.id_r]
    v = x[veh.id_v]
    quat = x[veh.id_quat]
    ω = x[veh.id_ω]

    if integrator.t >= integrator.p.MotorIgnitionTime[]
        traj = RocketTrajectoryParameters(r0=r, v0=v, q0=quat, ω0=ω, t0=integrator.t - integrator.p.MotorIgnitionTime[], MotorFired=true)
    else
        traj = RocketTrajectoryParameters(r0=r, v0=v, q0=quat, ω0=ω)
    end

    mdl = RocketProblem(veh, atmos, traj)
    tmpSolution = solveProblem(mdl)

    println(tmpSolution.status)

    println("\n", traj, "\n")

    if tmpSolution.status == "SCP_SOLVED"
        integrator.p.solution[] = tmpSolution
        
        if !traj.MotorFired
            integrator.p.MotorIgnitionTime[] = integrator.t + integrator.p.solution[].p[veh.id_tcoast]
        end
    end

    println(integrator.t + integrator.p.solution[].p[veh.id_tcoast])
end

mpccb = PeriodicCallback(mpc!, 0.25);

import TVC: motorTime
function motorTime(t, MotorIgnitionTime)
    return t - MotorIgnitionTime[]
end

@with_kw struct ODEParameters{R, S, V}
    veh::RocketParameters = RocketParameters()
    atmos::Atmosphere = Atmosphere()
    traj::RocketTrajectoryParameters = RocketTrajectoryParameters()

    Aero::Bool = false
    wind::V = zeros(3)

    solution::S

    MotorIgnitionTime::R = 0.0
    Control = zeros(3)
end

MotorIgnitionTime = Ref(solution.p[veh.id_tcoast])

tspan = (0., solution.p[veh.id_tcoast] + veh.BurnTime + 10.0)

test = Ref(solution)
p = ODEParameters(veh=veh, atmos=atmos, traj=traj, Aero=false, solution=test, MotorIgnitionTime=MotorIgnitionTime, Control=(x, p, t) -> control(x, p, t, servoΔt))

prob = ODEProblem(f!, x_0, tspan, p)

condition(x,t,integrator) = x[3] # when zero halt integration
cb = ContinuousCallback(condition, nothing, terminate!) # when going upwards do nothing

saved_values = SavedValues(Float64, Float64)
scb = SavingCallback((u, t, integrator) -> integrator.p.MotorIgnitionTime[], saved_values, save_everystep = true)

cbs = CallbackSet(cb, mpccb, scb);
# cbs = CallbackSet(cb, scb);

sol = DifferentialEquations.solve(prob, reltol=1e-8, abstol=1e-8, callback=cbs)

println(sol.prob.p.MotorIgnitionTime[])
println(sol.t[end] - sol.prob.p.MotorIgnitionTime[])
Plots.plot(sol, vars=1:3)
Plots.plot(sol, vars=4:6)
Plots.plot(sol, vars=7:10)

Plots.plot(saved_values.t, saved_values.saveval) 