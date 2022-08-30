using TVC, DifferentialEquations, SCPToolbox, Plots, Revise, LinearAlgebra, Parameters

#   Specify Parameters
#   ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡

veh = RocketParameters()
atmos = Atmosphere()
traj = RocketTrajectoryParameters();

mdl = RocketProblem(veh, atmos, traj)

# servoΔt = 0.02 # Servo step rate
x_0 = [traj.r0; traj.v0; traj.q0; traj.ω0][[veh.id_r; veh.id_v; veh.id_quat; veh.id_ω]];

#   Solve Optimal Trajectory
#   ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡

solution = solveProblem(mdl);

#   Specify Controller
#   ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡

# Fix mpc! before using.
# function control(x, p, t, servoΔt)
#     veh = p.veh
#     tₘ = motorTime(t, p.MotorIgnitionTime)
#     sol = p.solution[]
#     t0 = p.t0[]

#     td = tₘ - tₘ % servoΔt

#     if t0 ≤ tₘ ≤ veh.BurnTime
#         time = (td - t0) / (veh.BurnTime - t0)
#         desired_tvc = sample(sol.xc, time)[veh.id_T]
#         desired_roll = sample(sol.uc, time)[veh.id_roll]
#     else
#         desired_tvc = zeros(3)
#         desired_roll = 0.
#     end

#     return Actuator(x, p, t, desired_tvc, desired_roll)
# end

function control(x, p, t)
    veh = p.veh
    tₘ = motorTime(t, p.MotorIgnitionTime)
    sol = p.solution[]
    t0 = p.t0[]

    if t0 ≤ tₘ ≤ veh.BurnTime
        time = (tₘ - t0) / (veh.BurnTime - t0)
        desired_tvc = normalize(sample(sol.xc, time)[veh.id_T])
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

#   MPC Controller
#   ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡

function mpc!(integrator) # what if motor is spent, is this dealt with properly?
    veh = integrator.p.veh
    atmos = integrator.p.atmos

    x = integrator.u

    r = x[veh.id_r]
    v = x[veh.id_v]
    quat = x[veh.id_quat]
    ω = x[veh.id_ω]

    # Only valid for continuous control
    tₘ = motorTime(integrator.t, integrator.p.MotorIgnitionTime[])
    sol = p.solution[]
    t0 = p.t0[]

    if integrator.t >= integrator.p.MotorIgnitionTime[]
        time = (tₘ - t0) / (veh.BurnTime - t0)
        traj = RocketTrajectoryParameters(r0=r, v0=v, q0=quat, ω0=ω, T0=sample(sol.xc, time)[veh.id_T], Ṫ0=sample(sol.xc, time)[veh.id_Ṫ], t0=tₘ, MotorFired=true, PreviousTrajectoryState=sol.xc, PreviousTrajectoryInput=sol.uc, PreviousTrajectoryCurrentTime=time, UsePreviousTrajectory=true)
    else
        traj = RocketTrajectoryParameters(r0=r, v0=v, q0=quat, ω0=ω, PreviousTrajectoryState=sol.xc, PreviousTrajectoryInput=sol.uc, PreviousTrajectoryP=integrator.p.MotorIgnitionTime[] -integrator.t, UsePreviousTrajectory=true)
    end

    mdl = RocketProblem(veh, atmos, traj)
    tmpSolution = solveProblem(mdl)

    println(tmpSolution.status)

    println("\n", traj, "\n")

    if tmpSolution.status == "SCP_SOLVED"
        integrator.p.solution[] = tmpSolution  
        
        if !traj.MotorFired
            integrator.p.MotorIgnitionTime[] = integrator.t + integrator.p.solution[].p[veh.id_tcoast] # update motor ignition time as long as motor hasn't been fired
            # integrator.p.t0[] = 0. # unnecessary
        else
            integrator.p.t0[] = integrator.t - integrator.p.MotorIgnitionTime[] # t0 is 0 until motor has been fired, so once motor fired start updating.
        end
    end

    println(integrator.t)
    println(integrator.t + integrator.p.solution[].p[veh.id_tcoast])
    println()
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

    # mdl::RocketProblem = RcoketProblem()
    
    Aero::Bool = false
    wind::V = zeros(3)

    ground::Bool = false

    solution::S

    MotorIgnitionTime::R = 0.0
    t0::Base.RefValue{Float64} = Ref(0.0) # Seems to store time since motor started firing when guidance is solved.
    # confusing that traj and odeparamaters both have t0
    
    Control = (x, p, t) -> (force=[0; 0; veh.Thrust(t)], torque=zeros(3))
end

MotorIgnitionTime = Ref(solution.p[veh.id_tcoast])

tspan = (traj.t0, traj.t0 + solution.p[veh.id_tcoast] + veh.BurnTime + 10.0)

p = ODEParameters(veh=veh, atmos=atmos, traj=traj, Aero=true, solution=Ref(solution), MotorIgnitionTime=MotorIgnitionTime, Control=(x, p, t) -> control(x, p, t))

prob = ODEProblem(f!, x_0, tspan, p)

condition(x,t,integrator) = x[3] # when zero halt integration
cb = ContinuousCallback(condition, nothing, terminate!) # when going upwards do nothing

# saved_values = SavedValues(Float64, Float64)
# scb = SavingCallback((u, t, integrator) -> integrator.p.MotorIgnitionTime[], saved_values, save_everystep = true)

cbs = CallbackSet(cb, mpccb)#, scb);
# cbs = CallbackSet(cb, scb);

sol = DifferentialEquations.solve(prob, callback=cbs)
# sol = DifferentialEquations.solve(prob, reltol=1e-8, abstol=1e-8, callback=cbs)

println(sol.u[end][mdl.veh.id_v])

println(sol.prob.p.MotorIgnitionTime[])
println(sol.t[end] - sol.prob.p.MotorIgnitionTime[])

Plots.plot(sol, idxs=veh.id_r)
Plots.plot!(sol.prob.p.MotorIgnitionTime .+ sol.prob.p.traj.t0 .+ (mdl.veh.BurnTime - sol.prob.p.traj.t0) * solution.td, solution.xd[mdl.veh.id_r, :]')

Plots.plot(sol, idxs=4:6)
Plots.plot(sol, idxs=7:10)

# Plots.plot(saved_values.t, saved_values.saveval)