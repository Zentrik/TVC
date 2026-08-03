#=  Diagnostics for "the guidance problem often isn't solvable".

    See docs/mpc-feasibility.md for what these numbers mean. Briefly:

      * `statusSweep` shows that solvability is not monotone in the initial
        condition, which is the signature of a numerical breakdown in the
        conic solver rather than of an infeasible problem.
      * `lastSafeIterate` shows that when PTR reports SCP_FAILED it has almost
        always already produced a usable trajectory on an earlier iteration.
      * `verticalAuthority` shows how little of the terminal altitude
        constraint the guidance can still influence once the motor is lit.
      * `min‖T‖` in the summary lines is the tell for whether `‖T‖ ≤ 1` came
        out tight. Anything below 1 is a throttle setting a solid motor cannot
        produce, i.e. the relaxation was lossy and the plan is not flyable.
=#

using TVC, SCPToolbox, LinearAlgebra, Printf, ECOS
import JuMP: MOI

veh = RocketParameters()
atmos = Atmosphere()

"""
    solveWithHistory(mdl; kwargs...)

Same problem `solveProblem` builds, but returns `(sol, hist)` so the individual
PTR iterations can be inspected. `run.jl` throws the history away.
"""
function solveWithHistory(mdl; iter_max=50, ε_abs=1e-5, ε_rel=1e-3, solver=ECOS,
                          options=Dict{String, Any}("verbose" => 0))
    pbm = TrajectoryProblem(mdl)
    TVC.Guidance.define_problem!(pbm, :ptr)

    N = max(floor(Int, (mdl.veh.BurnTime - mdl.traj.t0) / 0.4) + 2, 5)
    pars = PTR.Parameters(N, 100, iter_max, FOH, 5e3, 1e-2, ε_abs, ε_rel, 1e-2,
                          Inf, Inf, solver, options)

    return PTR.solve(PTR.create(pars, pbm))
end

subproblemSolved(sol) = sol.status == MOI.OPTIMAL || sol.status == MOI.ALMOST_OPTIMAL

"""
    lastSafeIterate(hist)

Index of the last PTR iteration whose conic subproblem actually solved, or
`nothing`. `SCPToolbox` only looks at the last iteration when it decides
between `SCP_SOLVED` and `SCP_FAILED`, so a single bad final subproblem
discards every good iterate before it. For MPC that is the wrong trade: keep
this one instead.
"""
function lastSafeIterate(hist)
    solutions = [subproblem.sol for subproblem in hist.subproblems]
    return findlast(subproblemSolved, solutions)
end

function summarise(label, mdl; kwargs...)
    try
        sol, hist = solveWithHistory(mdl; kwargs...)
        solutions = [subproblem.sol for subproblem in hist.subproblems]
        k = lastSafeIterate(hist)

        thrustMagnitude = isnothing(k) ? NaN :
            minimum(norm(solutions[k].xd[mdl.veh.id_T, j])
                    for j = 1:size(solutions[k].xd, 2))

        @printf("  %-28s %-26s iterations = %2d, min‖T‖ = %.4f, last usable iterate = %s\n",
                label, sol.status, length(solutions), thrustMagnitude,
                isnothing(k) ? "none" :
                @sprintf("%d (J = %.3e, max|vd| = %.1e, max|vbc| = %.1e)",
                         k, solutions[k].J, maximum(abs, solutions[k].vd),
                         maximum(abs, solutions[k].vtc)))

        return sol.status
    catch e # a subproblem the conic solver could not solve is re-discretised
            # before its status is checked, which can throw SingularException
        @printf("  %-28s THREW %s\n", label, first(split(sprint(showerror, e), '\n')))
        return "THREW"
    end
end

#   Solvability is not monotone in the initial condition
#   ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡

function statusSweep(; heights=16.0:1.0:36.0, vehicle=veh, kwargs...)
    statuses = String[]

    for height in heights
        traj = RocketTrajectoryParameters(r0=[20.0, -4.0, height])
        push!(statuses, summarise(@sprintf("h0 = %.0f m", height),
                                  RocketProblem(vehicle, atmos, traj); kwargs...))
    end

    @printf("  ==> %d/%d solved\n", count(==("SCP_SOLVED"), statuses), length(statuses))

    return statuses
end

println("Ignition altitude sweep, default PTR settings")
statusSweep()

println("\nSame sweep, stopping one iteration earlier (ε_rel = 1e-2)")
statusSweep(ε_rel=1e-2) # no effect, the convergence test cannot see the bad
                        # subproblem coming

println("\nSame sweep with the old formulation: ‖T‖ ≤ 1, touchdown at burnout")
statusSweep(vehicle=RocketParameters(Throttleable=true, FixedLandingTime=true))

# `] add Clarabel` and uncomment: this solves every case ECOS fails on, though
# it needs more SCP iterations to get there.
# using Clarabel
# println("\nSame sweep, Clarabel instead of ECOS")
# statusSweep(solver=Clarabel, options=Dict{String, Any}("verbose" => false))

#   Re-solving part way through the burn
#   ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡

function midBurnSweep(; motorTimes=[0.5, 1.5, 2.5, 3.0],
                      altitudeErrors=[0.0, -0.5, 0.5, -2.0, 2.0])
    nominal, = solveWithHistory(RocketProblem(veh, atmos, RocketTrajectoryParameters()))

    for tₘ in motorTimes
        x = sample(nominal.xc, tₘ / veh.BurnTime)

        for Δh in altitudeErrors
            traj = RocketTrajectoryParameters(r0=x[veh.id_r] + [0; 0; Δh],
                                              v0=x[veh.id_v], q0=x[veh.id_quat],
                                              ω0=x[veh.id_ω], T0=x[veh.id_T],
                                              Ṫ0=x[veh.id_Ṫ], t0=tₘ, MotorFired=true)

            summarise(@sprintf("t0 = %.2f s, Δh = %+.1f m", tₘ, Δh),
                      RocketProblem(veh, atmos, traj))
        end
    end
end

println("\nRestarting from a point on the nominal trajectory, with an altitude error")
midBurnSweep()

#   How much of the terminal altitude constraint can still be influenced
#   ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡

"""
    verticalAuthority(t0)

Integrate the vertical channel from motor time `t0` to burnout at full throttle
and at the lowest throttle the `‖Ṫ‖ ≤ 5 deg/s` and `‖T̈‖ ≤ 10 deg/s²` limits
allow, starting from `‖T‖ = 1`. The gap between the two is the entire altitude
and velocity error the guidance can still absorb, given that the terminal
condition pins touchdown to burnout.
"""
function verticalAuthority(t0; Δt=1e-4)
    rate, acceleration = deg2rad(5), deg2rad(10)
    rampTime = rate / acceleration
    minThrottle = t -> max(0.0, 1 - (t - t0 <= rampTime ?
                                     acceleration * (t - t0)^2 / 2 :
                                     rate * (t - t0 - rampTime / 2)))

    function propagate(throttle)
        n = max(round(Int, (veh.BurnTime - t0) / Δt), 1)
        h = (veh.BurnTime - t0) / n
        height, velocity = 0.0, 0.0

        for i = 1:n
            t = t0 + (i - 0.5) * h
            a = throttle(t) * veh.Acceleration(t) - 9.80655
            height += velocity * h + a * h^2 / 2
            velocity += a * h
        end

        return height, velocity
    end

    return propagate(t -> 1.0) .- propagate(minThrottle)
end

println("\nRemaining vertical authority at burnout, as a function of when we re-solve")
println("  burn left | altitude window | velocity window")

for t0 = 0.0:0.5:3.0
    Δh, Δv = verticalAuthority(t0)
    @printf("   %5.2f s  | %10.3f m   | %8.3f m/s\n", veh.BurnTime - t0, Δh, Δv)
end
