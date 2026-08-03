# TODO: Make quaternion discretisation more accurate. Without constraint on quaternion norm, the norm of the quaternion strays from 1, I think this is due to the Jacobians of dynamics being evaluated about a reference trajectory that has non unit quaternions.

using SCPToolbox
# import SCPToolbox.Parser.@perturb_fix

using LinearAlgebra
# using ECOS

using ..Utils
import ..Utils: rotate, skew
using ..Guidance

using ForwardDiff
using JuMP

export define_problem!

# Per Successive Convexification for Mars 6-DoF Powered Descent Landing Guidance, 2017. Set control to second derivative of thrust vector
# This allows me to add a constraint on the tvc gimbal rate and it adds more degrees of freedom on the control/ allows for more complex control for a given N (no. of discretisation steps) improving cost.

function define_problem!(pbm::TrajectoryProblem, algo::Symbol)::Nothing
    set_dims!(pbm)
    set_scale!(pbm)
    set_cost!(pbm)
    set_dynamics!(pbm)
    set_integration_action(pbm)
    set_convex_constraints!(pbm)
    set_nonconvex_constraints!(pbm, algo)
    set_bcs!(pbm)

    set_guess!(pbm)

    return nothing
end

function set_dims!(pbm::TrajectoryProblem)::Nothing

    problem_set_dims!(pbm, 19, 4, 2) # parameters are [t_coast; t_land]

    return nothing
end

function set_scale!(pbm::TrajectoryProblem)::Nothing #VERY IMPORTANT
    advise! = problem_advise_scale!

    # States
    advise!(pbm, :state, 1, (-100.0, 100.0))
    advise!(pbm, :state, 2, (-100.0, 100.0))
    advise!(pbm, :state, 3, (0.0, 100.0))
    advise!(pbm, :state, 4, (-50.0, 50.0))
    advise!(pbm, :state, 5, (-50.0, 50.0))
    advise!(pbm, :state, 6, (-50.0, 50.0))

    advise!(pbm, :state, 7, (-1.0, 1.0))
    advise!(pbm, :state, 8, (-1.0, 1.0))
    advise!(pbm, :state, 9, (-1.0, 1.0))
    advise!(pbm, :state, 10, (-1.0, 1.0))
    advise!(pbm, :state, 11, (-10.0, 10.0))
    advise!(pbm, :state, 12, (-10.0, 10.0))
    advise!(pbm, :state, 13, (-10.0, 10.0)) # was (-10.0, 0.0), which is a typo: ω_z is not sign definite

    advise!(pbm, :state, 14, (-1.0, 1.0))
    advise!(pbm, :state, 15, (-1.0, 1.0))
    advise!(pbm, :state, 16, (-1.0, 1.0))
    advise!(pbm, :state, 17, (-float(π), float(π)))
    advise!(pbm, :state, 18, (-float(π), float(π)))
    advise!(pbm, :state, 19, (-float(π), float(π)))

    # Inputs
    advise!(pbm, :input, 1, (-deg2rad(10), deg2rad(10)))
    advise!(pbm, :input, 2, (-deg2rad(10), deg2rad(10)))
    advise!(pbm, :input, 3, (-deg2rad(10), deg2rad(10)))
    advise!(pbm, :input, 4, (-1.0, 1.0))

    # Parameters
    advise!(pbm, :parameter, 1, (0.0, 10.0))
    advise!(pbm, :parameter, 2, (0.0, pbm.mdl.veh.BurnTime))

    return nothing
end

function straightline_interpolate(v0, vf, N::Int)

    # Initialize
    nv = length(v0)
    v = zeros(nv, N)

    for k = 1:N
        mix = (k - 1) / (N - 1)

        v[:, k] = v0 + mix * (vf - v0)
    end

    return v
end

function set_guess!(pbm::TrajectoryProblem)::Nothing
    problem_set_guess!(pbm, (N, pbm) -> begin
        veh = pbm.mdl.veh
        traj = pbm.mdl.traj
        atmos = pbm.mdl.atmos

        if traj.UsePreviousTrajectory
            p = [traj.PreviousTrajectoryP; traj.PreviousTrajectoryTLand]

            SampleTimes = collect(range(traj.PreviousTrajectoryCurrentTime, 1, N))
            x = mapreduce(t -> sample(traj.PreviousTrajectoryState, t), hcat, SampleTimes)
            u = mapreduce(t -> sample(traj.PreviousTrajectoryInput, t), hcat, SampleTimes)
        else
            p = [0.0; veh.BurnTime] # ignite immediately, land at burnout.

            # motorTimeRemaining = veh.BurnTime - traj.t0 # how much motor time remaining

            # dt = motorTimeRemaining / (N - 1) # for convex constraints

            x = zeros(pbm.nx, N)
            u = zeros(pbm.nu, N)

            x[veh.id_r, :] = straightline_interpolate(traj.r0, [traj.r0[1:2]; traj.rN[3]], N)
            x[veh.id_v, :] = straightline_interpolate(traj.v0, traj.vN, N)
            x[veh.id_ω, :] = straightline_interpolate(traj.ω0, traj.ωN, N)

            x[veh.id_T, :] = straightline_interpolate(traj.T0, traj.TN, N)  
            x[veh.id_Ṫ, :] = straightline_interpolate(traj.Ṫ0, traj.ṪN, N) 
            
            u = straightline_interpolate([0; 0; 0; 0], [0; 0; 0; 0], N)
            
            # Roll is left at whatever q0 has, only the tilt is interpolated out.
            # This has to be normalised: [q0[1]; qN[2:3]; q0[4]] has norm < 1
            # whenever the rocket is tilted, and a non unit reference quaternion
            # makes the linearised dynamics inconsistent (see the TODO at the
            # top of this file).
            quatN = normalize([traj.q0[1]; traj.qN[2:3]; traj.q0[4]])

            for k = 1:N
                mix = (k - 1) / (N - 1)

                x[veh.id_quat, k] = slerp_quat(traj.q0, quatN, mix)
            end
        end

        return x, u, p
    end)

    return nothing
end
    
function set_cost!(pbm::TrajectoryProblem)::Nothing
    problem_set_terminal_cost!(
        pbm, (x, p, pbm) -> dot(x[pbm.mdl.veh.id_v] - pbm.mdl.traj.vN, x[pbm.mdl.veh.id_v] - pbm.mdl.traj.vN)
        # 0 # use for feasibility testing
    )

    return nothing
end

function set_dynamics!(pbm::TrajectoryProblem)::Nothing
    function f_t(x, u, t, pbm)
        veh = pbm.mdl.veh
        atmos = pbm.mdl.atmos

        r = x[veh.id_r]
        v = x[veh.id_v]
        quat = x[veh.id_quat]
        ω = x[veh.id_ω]
        T = x[veh.id_T]
        Ṫ = x[veh.id_Ṫ]

        return [v; atmos.g(r[3]) + rotate(quat, T) * veh.Acceleration(t); 1/2 * quatL(quat) * [0; ω]; veh.InertiaTensor \ (cross(veh.MomentArm(t), T * veh.Thrust(t)) + [0; 0; u[veh.id_roll]] - cross(ω, veh.InertiaTensor * ω)); Ṫ; u[veh.id_T̈]]
    end

    function A_t(x, u, t, veh)
        r = x[veh.id_r]
        v = x[veh.id_v]
        quat = x[veh.id_quat]
        ω = x[veh.id_ω]
        T = x[veh.id_T]
        Ṫ = x[veh.id_Ṫ]

        A = zeros(pbm.nx, pbm.nx)
        A[veh.id_r, veh.id_v] = I(3)
        A[veh.id_v, veh.id_quat[1]] = veh.Acceleration(t) * 2 * (quat[1] * T + quat[2:4] × T)
        A[veh.id_v, veh.id_quat[2:4]] = veh.Acceleration(t) * 2 * (quat[2:4]' * T * I(3) + quat[2:4] * T' - T * quat[2:4]' - quat[1] * skew(T))
        A[veh.id_v, veh.id_T] = veh.Acceleration(t) * to_matrix(quat)
        A[veh.id_quat, veh.id_quat] = 1/2 * quatR([0; ω])
        A[veh.id_quat, veh.id_ω] = 1/2 * quatL(quat)[:, 2:4]
        A[veh.id_ω, veh.id_ω] = - veh.InertiaTensor \ (skew(ω) * veh.InertiaTensor - skew(veh.InertiaTensor * ω))
        A[veh.id_ω, veh.id_T] = veh.InertiaTensor \ skew(veh.MomentArm(t)) * veh.Thrust(t)
        A[veh.id_T, veh.id_Ṫ] = I(3)

        return A
    end

    function B_t(x, u, t, veh)
        B = zeros(pbm.nx, pbm.nu)
        B[veh.id_ω, veh.id_roll] = inv(veh.InertiaTensor)[:, 3]
        B[veh.id_Ṫ, veh.id_T̈] = I(3)

        return B
    end

    # ∂f/∂t. `Thrust`, `Mass` and `CG` are interpolated tables, so this is done
    # numerically rather than by differentiating through them. It is only needed
    # for df/dp: stretching the horizon slides every node along the thrust curve
    # as well as changing dt/dτ.
    function ḟ_t(x, u, t, pbm, h=1e-5)
        return (f_t(x, u, t + h, pbm) - f_t(x, u, t - h, pbm)) / (2 * h)
    end

    # Dynamics, in normalised time τ ∈ [0, 1]. The horizon is t_land - t0, where
    # t_land is a decision variable, so dt/dτ = t_land - t0 and every derivative
    # below picks up that factor.
    problem_set_dynamics!(
        pbm,
        # f
        (t, k, x, u, p, pbm) -> f_t(x, u, motorTime(t, p, pbm.mdl), pbm) * horizon(p, pbm.mdl),
        # df/dx
        (t, k, x, u, p, pbm) -> A_t(x, u, motorTime(t, p, pbm.mdl), pbm.mdl.veh) * horizon(p, pbm.mdl),
        # df/du
        (t, k, x, u, p, pbm) -> B_t(x, u, motorTime(t, p, pbm.mdl), pbm.mdl.veh) * horizon(p, pbm.mdl),
        # df/dp
        (t, k, x, u, p, pbm) -> begin
            veh = pbm.mdl.veh
            tₘ = motorTime(t, p, pbm.mdl)

            F = zeros(pbm.nx, pbm.np)
            # d/dt_land of f(x, u, t0 + τ (t_land - t0)) * (t_land - t0)
            F[:, veh.id_tland] = f_t(x, u, tₘ, pbm) +
                                 horizon(p, pbm.mdl) * t * ḟ_t(x, u, tₘ, pbm)

            F
        end
    )

    return nothing
end

function set_integration_action(pbm::TrajectoryProblem)::Nothing

    # Quaternion re-normalization on numerical integration step
    problem_set_integration_action!(
        pbm, pbm.mdl.veh.id_quat,
        (q, pbm) -> begin
            if norm(q) ≈ 0
                print(q)
                qn = q
            else
                qn = q / norm(q)
            end
            
            return qn
    end)

    return nothing
end

function set_bcs!(pbm::TrajectoryProblem)::Nothing
    # Boundary conditions
    
    if false #pbm.mdl.traj.MotorFired # If motor has been fired, coast time is 0
        problem_set_bc!(
            pbm, :ic, # Initial condition
            (x, p, pbm) -> begin
                veh = pbm.mdl.veh
                traj = pbm.mdl.traj

                x0 = zeros(pbm.nx)
                x0[veh.id_r] .= traj.r0
                x0[veh.id_v] .= traj.v0
                x0[veh.id_quat] .= traj.q0
                x0[veh.id_ω] = traj.ω0
                x0[veh.id_T] = traj.T0
                x0[veh.id_Ṫ] = traj.Ṫ0

                return x - x0
            end,
            (x, p, pbm) -> I(pbm.nx), # Jacobian wrt x
            (x, p, pbm) -> zeros(pbm.nx, pbm.np) # Jacobian wrt p
        )
    else
        problem_set_bc!(
            pbm, :ic, # Initial condition
            (x, p, pbm) -> begin
                veh = pbm.mdl.veh
                traj = pbm.mdl.traj
                atmos = pbm.mdl.atmos

                x0 = zeros(pbm.nx)
                x0[veh.id_r] = traj.r0 + traj.v0 * p[veh.id_tcoast] + atmos.g(traj.r0[3]) * p[veh.id_tcoast]^2/2 # is there a better approximation? can we just calculate the exact result. (rn g is constant, so its fine)
                x0[veh.id_v] = traj.v0 + atmos.g(traj.r0[3]) * p[veh.id_tcoast]
                x0[veh.id_quat] = quatL(traj.q0) * wexp(traj.ω0 * p[veh.id_tcoast])
                x0[veh.id_ω] = traj.ω0 # not true is it?
                x0[veh.id_T] = traj.T0 + traj.Ṫ0 * p[veh.id_tcoast] # this should be fine
                x0[veh.id_Ṫ] = traj.Ṫ0 # this should be true

                if  p[veh.id_tcoast] > 1e-3 && norm(veh.InertiaTensor \ (- cross(traj.ω0, veh.InertiaTensor * traj.ω0))) > 1e-5 #abs(traj.ω0[3]) > 1e-4
                    println("Warning: Initial constraint probably false, angular acceleration in free fall is non zero") 
                end

                return x - x0
            end,
            (x, p, pbm) -> I(pbm.nx), # Jacobian wrt x
            (x, p, pbm) -> begin # Jacobian wrt p 
                veh = pbm.mdl.veh
                traj = pbm.mdl.traj
                atmos = pbm.mdl.atmos

                J = zeros(pbm.nx, pbm.np)
                J[veh.id_r, veh.id_tcoast] = traj.v0 + atmos.g(traj.r0[3]) * p[veh.id_tcoast]
                J[veh.id_v, veh.id_tcoast] = atmos.g(traj.r0[3])
                J[veh.id_quat, veh.id_tcoast] = ForwardDiff.derivative(t -> quatL(traj.q0) * wexp(traj.ω0 * t), p[veh.id_tcoast])
                # J[veh.id_q, veh.id_tcoast] = quatL(traj.q0) * [-sin(norm(traj.ω0 * p[veh.id_tcoast]) / 2); traj.ω0 / norm(traj.ω0) * cos(norm(traj.ω0 * p[veh.id_tcoast]) / 2)] * norm(traj.ω0) / 2
                # will fail if norm(traj.ω0) = 0
                J[veh.id_ω, veh.id_tcoast] = zeros(3)
                J[veh.id_T, veh.id_tcoast] = traj.Ṫ0
                J[veh.id_Ṫ, veh.id_tcoast] = zeros(3)

                # if norm(ForwardDiff.derivative(t -> quatL(traj.q0) * wexp(traj.ω0 * t), p[veh.id_tcoast]) - J[veh.id_q, veh.id_tcoast]) > 1e-5
                #     println("Warning: Intial costraint derivative wrong.") 
                # end 

                return -J
            end,
        )
    end

    problem_set_bc!(
            pbm, :tc, # Terminal condition
            (x, p, pbm) -> begin
                veh = pbm.mdl.veh
                traj = pbm.mdl.traj

                xf = zeros(6)
                xf[1] = traj.rN[3]
                xf[2:3] = traj.qN[2:3]
                xf[4:6] = traj.ωN

                return x[vcat(veh.id_r[3], veh.id_quat[2:3], veh.id_ω)] - xf
            end,
            (x, p, pbm) -> begin # Jacobian wrt x 
                veh = pbm.mdl.veh
                traj = pbm.mdl.traj

                J = zeros(6, pbm.nx)
                J[1, veh.id_r[3]] = 1
                J[2:3, veh.id_quat[2:3]] = I(2)
                J[4:6, veh.id_ω] = I(3)

                return J
            end,
            (x, p, pbm) -> zeros(6, pbm.np), # Jacobian wrt p
        )

        # problem_set_bc!(
        #     pbm, :tc, # Terminal condition
        #     (x, p, pbm) -> begin
        #         veh = pbm.mdl.veh
        #         traj = pbm.mdl.traj

        #         xf = zeros(9)
        #         xf[1] = traj.rN[3]
        #         xf[2:3] = traj.qN[2:3]
        #         xf[4:6] = traj.ωN
        #         xf[7:9] = traj.vN

        #         return x[vcat(veh.id_r[3], veh.id_quat[2:3], veh.id_ω, veh.id_v)] - xf
        #     end,
        #     (x, p, pbm) -> begin # Jacobian wrt x 
        #         veh = pbm.mdl.veh
        #         traj = pbm.mdl.traj

        #         J = zeros(9, pbm.nx)
        #         J[1, veh.id_r[3]] = 1
        #         J[2:3, veh.id_quat[2:3]] = I(2)
        #         J[4:6, veh.id_ω] = I(3)
        #         J[7:9, veh.id_v] = I(3)

        #         return J
        #     end,
        #     (x, p, pbm) -> zeros(9, pbm.np), # Jacobian wrt p
        # )
end

"""Shortest powered horizon the guidance will plan over, in seconds."""
const MinimumHorizon = 0.05

"""
    fixParameter!(ocp, pbm, i, value)

Pin element `i` of the problem's parameter vector to the physical value `value`.

Two things to be careful of. The parameter block is created as a vector, so JuMP
names its elements `p[1]`, `p[2]`, ... — asking for `"p"` gets you `nothing` as
soon as there is more than one of them. And the JuMP variable is the *scaled*
parameter, `p = S p̂ + c` with `S`, `c` taken from the bounds given to
`problem_advise_scale!`, so the value has to be scaled to match.
"""
function fixParameter!(ocp, pbm::TrajectoryProblem, i::Integer, value::Real)
    lower, upper = pbm.prg[i]

    fix(variable_by_name(jump_model(ocp), "p[$i]"), (value - lower) / (upper - lower))

    return nothing
end

function set_convex_constraints!(pbm::TrajectoryProblem)::Nothing
    # Convex State Constraints
    problem_set_X!(
        pbm, (t, k, x, p, pbm, ocp) -> begin
            veh = pbm.mdl.veh
            traj = pbm.mdl.traj

            @add_constraint(
                ocp, NONPOS, "height >= 0", (x[veh.id_r[3]],), begin
                local height = arg[1]
                - height
                end)

            if traj.MotorFired
                # @add_constraint(
                #     ocp, ZERO, "t_coast == 0", (p[veh.id_tcoast],), begin
                #         local t_coast = arg[1]
                #         t_coast
                #     end)
                # `p` is a vector, so JuMP names its elements p[1], p[2], ...
                fixParameter!(ocp, pbm, veh.id_tcoast, 0.)
                # Probably doesn't matter, Gurobi seems to be able to equate the two above in its presolve, ECOS doesn't but it returns p on the order of 1e-8 or below with first constraint.

                # @perturb_fix p[veh.id_tcoast] # fix to initial guess which is 0? # doesn't seem to work well
            else
                @add_constraint(
                    ocp, NONPOS, "t_coast >= expected ignition time", (p[veh.id_tcoast],), begin
                        local t_coast = arg[1]
                        traj.ExpectedIgnitionTime - t_coast
                    end)
            end

            if veh.FixedLandingTime
                fixParameter!(ocp, pbm, veh.id_tland, veh.BurnTime)
            else
                # t0 + MinimumHorizon ≤ t_land ≤ BurnTime. The motor cannot be
                # relit and `veh.Thrust` is 0 past BurnTime, so there is nothing
                # to plan with beyond it; the lower bound just keeps the horizon
                # (and hence the scaled dynamics) away from zero.
                @add_constraint(
                    ocp, NONPOS, "t_land <= BurnTime", (p[veh.id_tland],), begin
                        local t_land = arg[1]
                        t_land - veh.BurnTime
                    end)

                @add_constraint(
                    ocp, NONPOS, "t_land >= t0 + minimum horizon", (p[veh.id_tland],), begin
                        local t_land = arg[1]
                        # `min` so this can never contradict the bound above,
                        # however little burn time is left.
                        min(traj.t0 + MinimumHorizon, veh.BurnTime) - t_land
                    end)
            end

            @add_constraint(
                ocp, SOC, "Thrust Magnitude <= Max", (x[veh.id_T],), begin # we have say thrust <= 1, as we want it normalised
                    local Thrust = arg[1]
                    [1; Thrust]
                end) # if the motor cannot throttle this is only half of ‖T‖ = 1, see set_nonconvex_constraints!

            @add_constraint(
                ocp, SOC, "Thrust Gimal angle <= delta_max", (x[veh.id_T],), begin
                    local Thrust = arg[1]
                    [Thrust[3] / cos(pi * 5/180); Thrust]
                end)
            
            @add_constraint( # shouldn't be necessary, but seems very very useful.
                ocp, SOC, "|quat| <= 1", (x[veh.id_quat],), begin
                    local quat = arg[1]
                    [1; quat]
                end)

            @add_constraint(
            ocp, SOC, "TVC angular velocity <= Max", (x[veh.id_Ṫ],), (u) -> [deg2rad(5); u]) # Angular velocity is r × v / ||r||², v = u, assume u is ⊥ r and ||r||² = 1, so angular velocity is v = u.

        # @add_constraint(
        #     ocp, SOC, "|w| < w_max", (x[11:13],), begin
        #         local w = arg[1]
        #         [pi / 2; w]
        #     end)
    end)

    # Convex Input Constraints
    problem_set_U!(
        pbm, (t, k, u, p, pbm, ocp) -> begin
            veh = pbm.mdl.veh
            traj = pbm.mdl.traj

            @add_constraint(
                ocp, L1, "|roll torque| < max", (u[4],), begin
                    local u = arg[1]
                    [0.1; u]
                end)

            @add_constraint(
                ocp, SOC, "TVC Acceleration <= Max", (u[veh.id_T̈],), begin
                    local u = arg[1]
                    [deg2rad(10); u] # Angular Acceleration is r × a / ||r||², a = u, assume u is ⊥ r and ||r||² = 1, so angular acceleration is a = u.
                end)
    end)
end

"""
    set_nonconvex_constraints!(pbm, algo)

A solid motor cannot throttle, so the thrust magnitude is not a control: ‖T‖ = 1
for the whole burn. `set_convex_constraints!` already imposes ‖T‖ ≤ 1 as a
second order cone; the other half, ‖T‖ ≥ 1, is nonconvex, so it is handed to the
SCP algorithm as `s(x) = 1 - T ⋅ T ≤ 0` and linearised about the reference
trajectory. Like the dynamics it is relaxed with a (penalised) virtual control,
so it can never make a subproblem infeasible.

Note that the "TVC angular velocity" and "TVC Acceleration" constraints only
mean what their names say when ‖T‖ = 1 — they bound ‖Ṫ‖ and ‖T̈‖, which are the
gimbal rate and angular acceleration only for a unit thrust vector.

Set `RocketParameters(Throttleable=true)` to drop this and let the optimiser
pick ‖T‖ ∈ [0, 1] instead, which is what this problem used to do.
"""
function set_nonconvex_constraints!(pbm::TrajectoryProblem, algo::Symbol)::Nothing
    if pbm.mdl.veh.Throttleable
        if algo == :scvx # SCvx wants an s even when there is nothing to enforce
            problem_set_s!(pbm, algo, (t, k, x, u, p, pbm) -> [0])
        end

        return nothing
    end

    problem_set_s!(
        pbm, algo,
        # s
        (t, k, x, u, p, pbm) -> begin
            local T = x[pbm.mdl.veh.id_T]

            [1 - dot(T, T)]
        end,
        # ds/dx
        (t, k, x, u, p, pbm) -> begin
            local C = zeros(1, pbm.nx)
            C[1, pbm.mdl.veh.id_T] = -2 * x[pbm.mdl.veh.id_T]

            C
        end,
        # ds/du
        (t, k, x, u, p, pbm) -> zeros(1, pbm.nu),
        # ds/dp
        (t, k, x, u, p, pbm) -> zeros(1, pbm.np),
    )

    return nothing
end