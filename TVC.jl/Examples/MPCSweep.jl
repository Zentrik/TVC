#=  Closed loop MPC, no DifferentialEquations.

    Plant  : TVC's own f!, fixed step RK4, Aero = true (the disturbance the
             guidance never sees).
    Control: thrust direction from the current plan, applied at full motor
             thrust — the vehicle cannot throttle.
    MPC    : solveProblem re-solved every 0.25 s from the measured state, with
             the same guards as Examples/MPC_Simulation.jl.
=#
using TVC, SCPToolbox, LinearAlgebra, Printf, Clarabel

const atmos = Atmosphere()
const MinimumReplanTime = 0.5
const RollRateGain = 10.0   # 1/s on the roll rate error
const MaxRollTorque = 0.1   # N m, matches the guidance problem's limit

mutable struct Plant
    veh::Any
    atmos::Any
    Aero::Bool
    wind::Vector{Float64}
    ground::Bool
    MotorIgnitionTime::Float64
    Control::Any
end
Plant(veh, ignition) = Plant(veh, atmos, true, zeros(3), false, ignition,
                             (x, p, t) -> (force = zeros(3), torque = zeros(3)))

function rk4(x, p, t, dt)
    k1 = zeros(13); k2 = zeros(13); k3 = zeros(13); k4 = zeros(13)
    f!(k1, x, p, t)
    f!(k2, x .+ dt / 2 .* k1, p, t + dt / 2)
    f!(k3, x .+ dt / 2 .* k2, p, t + dt / 2)
    f!(k4, x .+ dt .* k3, p, t + dt)
    y = x .+ dt / 6 .* (k1 .+ 2 .* k2 .+ 2 .* k3 .+ k4)
    y[7:10] ./= norm(y[7:10])
    return y
end

tiltAngle(q) = rad2deg(acos(clamp(TVC.Utils.to_matrix(q)[3, 3], -1, 1)))

function fly(label, traj0; veh = RocketParameters(), mpcΔt = 0.25, dt = 0.002)
    plan = solveProblem(RocketProblem(veh, atmos, traj0))
    if startswith(plan.status, "SCP_FAILED")
        @printf("%-26s initial solve failed: %s\n", label, plan.status)
        flush(stdout); return
    end

    ignition = plan.p[veh.id_tcoast]
    t0Plan = 0.0                       # motor time the current plan starts from
    plant = Plant(veh, ignition)
    x = [traj0.r0; traj0.v0; traj0.q0; traj0.ω0]

    function control(x, p, t)
        tₘ = t - p.MotorIgnitionTime
        tLand = plan.p[veh.id_tland]
        if t0Plan <= tₘ <= tLand && tₘ >= 0
            τ = clamp((tₘ - t0Plan) / (tLand - t0Plan), 0.0, 1.0)
            dir = normalize(sample(plan.xc, τ)[veh.id_T])
            # Feedforward the plan's roll torque and close a loop on the roll
            # *rate error*. The plan's ω_z is what its gyroscopic steering needs
            # and it ends at zero, so tracking it keeps the effector while
            # actually arriving with the roll the plan intended. Driving ω_z to
            # zero instead, or bounding it in the guidance, removes the effector
            # and the vehicle tumbles — see docs/mpc-feasibility.md.
            roll = clamp(sample(plan.uc, τ)[veh.id_roll][1] +
                         RollRateGain * veh.InertiaTensor[3, 3] *
                         (sample(plan.xc, τ)[veh.id_ω[3]] - x[veh.id_ω[3]]),
                         -MaxRollTorque, MaxRollTorque)
        else
            dir = zeros(3); roll = 0.0
        end


        F = dir * veh.Thrust(tₘ)
        return (force = F, torque = veh.MomentArm(tₘ) × F + [0; 0; roll])
    end
    plant.Control = control

    t, nextMpc, solves, failures = 0.0, mpcΔt, 0, 0
    tMax = ignition + veh.BurnTime + veh.MaxBallisticTime + 5.0

    while t < tMax
        if t >= nextMpc - 1e-12
            nextMpc += mpcΔt
            tₘ = t - plant.MotorIgnitionTime

            if tₘ <= veh.BurnTime - MinimumReplanTime
                solves += 1
                r, v, q, ω = x[1:3], x[4:6], x[7:10], x[11:13]
                tLand = plan.p[veh.id_tland]
                if tₘ >= 0
                    τ = clamp((tₘ - t0Plan) / (tLand - t0Plan), 0.0, 1.0)
                    xs = sample(plan.xc, τ)
                    tr = RocketTrajectoryParameters(r0 = r, v0 = v, q0 = q, ω0 = ω,
                        T0 = xs[veh.id_T], Ṫ0 = xs[veh.id_Ṫ], t0 = tₘ, MotorFired = true,
                        PreviousTrajectoryState = plan.xc, PreviousTrajectoryInput = plan.uc,
                        PreviousTrajectoryCurrentTime = τ, PreviousTrajectoryTLand = tLand,
                        UsePreviousTrajectory = true)
                else
                    tr = RocketTrajectoryParameters(r0 = r, v0 = v, q0 = q, ω0 = ω,
                        PreviousTrajectoryState = plan.xc, PreviousTrajectoryInput = plan.uc,
                        PreviousTrajectoryP = plant.MotorIgnitionTime - t,
                        PreviousTrajectoryTLand = tLand, UsePreviousTrajectory = true)
                end

                new = try
                    solveProblem(RocketProblem(veh, atmos, tr))
                catch
                    nothing
                end

                if !isnothing(new) && new.status == "SCP_SOLVED" && all(isfinite, new.xd)
                    plan = new
                    if tₘ < 0
                        plant.MotorIgnitionTime = t + new.p[veh.id_tcoast]
                    else
                        t0Plan = tₘ
                    end
                else
                    failures += 1
                end
            end
        end

        xNext = rk4(x, plant, t, dt)
        if xNext[3] <= 0 < x[3]                      # ground crossing
            α = x[3] / (x[3] - xNext[3])
            x = x .+ α .* (xNext .- x)
            t += α * dt
            break
        end
        x = xNext
        t += dt
    end

    tₘ = t - plant.MotorIgnitionTime
    v = x[4:6]
    @printf("%-26s |v|=%5.2f (vz=%6.2f, vxy=%5.2f) tilt=%5.1f° |ω|=%5.2f  landed t_m=%5.2f (burn %.2f)  r_xy=[%6.2f,%6.2f]  %d solves, %d rejected\n",
        label, norm(v), v[3], norm(v[1:2]), tiltAngle(x[7:10]), norm(x[11:13]),
        tₘ, veh.BurnTime, x[1], x[2], solves, failures)
    flush(stdout)
end

cases = [
    ("nominal",            RocketTrajectoryParameters()),
    ("5 m higher",         RocketTrajectoryParameters(r0 = [20.0, -4, 35])),
    ("5 m lower",          RocketTrajectoryParameters(r0 = [20.0, -4, 25])),
    ("faster horizontal",  RocketTrajectoryParameters(v0 = [6.0, -5, 0])),
    ("already descending", RocketTrajectoryParameters(v0 = [4.0, -3, -3])),
    ("released tilted 5°", RocketTrajectoryParameters(
        q0 = [cosd(2.5), sind(2.5), 0.0, 0.0], ω0 = [0.1, -0.1, 0.0])),
]

# One case per process (`julia MPCSweep.jl 3`) so a sweep can be run in
# parallel — a single flight is ~15 guidance solves and takes a few minutes.
if isempty(ARGS)
    for (label, traj) in cases
        fly(label, traj)
    end
else
    i = parse(Int, ARGS[1])
    fly(cases[i][1], cases[i][2])
end
