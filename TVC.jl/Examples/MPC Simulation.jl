using TVC, DifferentialEquations, SCPToolbox, Plots

veh = RocketParameters()
atmos = Atmosphere()
traj = RocketTrajectoryParameters();

mdl = RocketProblem(veh, atmos, traj)

servoΔt = 0.02 # Servo step rate
x_0 = [traj.r0; traj.v0; traj.q0; traj.ω0][[veh.id_r; veh.id_v; veh.id_quat; veh.id_ω]];

solution = solveProblem(mdl);

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
    integrator.p.solution[] = solveProblem(mdl)

    print(integrator.t + integrator.p.solution[].p[veh.id_tcoast])

    if !traj.MotorFired
        integrator.p.MotorIgnitionTime[] = integrator.t + integrator.p.solution[].p[veh.id_tcoast]
    end
end

mpccb = PeriodicCallback(mpc!, 1.);

import TVC: motorTime
function motorTime(t, MotorIgnitionTime)
    return t - MotorIgnitionTime[]
end

function control(t, Rocket, solution, servoΔt) # discretise TVC based on servo step rate
    td = t - t % servoΔt
    if 0 ≤ t ≤ Rocket.BurnTime
        return vcat(sample(solution[].xc, td / Rocket.BurnTime)[Rocket.id_T] * Rocket.Thrust(t), sample(solution[].uc, t / Rocket.BurnTime)[Rocket.id_roll])
    else
        return zeros(4)
    end
end

MotorIgnitionTime = Ref(solution.p[veh.id_tcoast])

tspan = (0., solution.p[veh.id_tcoast] + veh.BurnTime + 10.0)
test = Ref(solution)
p = (veh=veh, atmos=atmos, wind=randn(3) * 5, MotorIgnitionTime=MotorIgnitionTime, solution=test, Control=t -> control(motorTime(t, MotorIgnitionTime), veh, test, servoΔt))

prob = ODEProblem(f!, x_0, tspan, p)

condition(x,t,integrator) = x[3] # when zero halt integration
cb = ContinuousCallback(condition, nothing, terminate!) # when going upwards do nothing

saved_values = SavedValues(Float64, Float64)
scb = SavingCallback((u, t, integrator) -> integrator.p.MotorIgnitionTime[], saved_values, save_everystep = true)

cbs = CallbackSet(cb, mpccb, scb);