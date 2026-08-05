# function scvx()
#     N = floor(Int, 3.45 / 0.1) + 1
#     Nsub = 100
#     iter_max = 30
#     disc_method = FOH
#     λ = 5e3
#     ρ_0 = 0.0
#     ρ_1 = 0.1
#     ρ_2 = 0.7
#     β_sh = 2.0
#     β_gr = 2.0
#     η_init = 1.0
#     η_lb = 1e-8
#     η_ub = 10.0
#     ε_abs = 1e-8
#     ε_rel = 1e-5
#     feas_tol = 5e-3
#     q_tr = Inf
#     q_exit = Inf
#     solver = ECOS
#     solver_options = Dict("verbose"=>0, "maxit"=>1000)
#     pars = SCvx.Parameters(N, Nsub, iter_max, disc_method, λ, ρ_0, ρ_1, ρ_2, β_sh, β_gr,
#                         η_init, η_lb, η_ub, ε_abs, ε_rel, feas_tol, q_tr,
#                         q_exit, solver, solver_options)


#     scvx_pbm = SCvx.create(pars, pbm)
#     sol, history = SCvx.solve(scvx_pbm)
#     return sol
# end

export solveProblem, lastUsableSolution

"""
    gridSize(mdl)

Number of discretisation nodes. Sized for the longest horizon the problem is
allowed to choose — burn time plus the ballistic tail — because the grid is
uniform in normalised time.
"""
function gridSize(mdl)
    horizon = mdl.veh.BurnTime + (mdl.veh.FixedLandingTime ? 0. : mdl.veh.MaxBallisticTime) - mdl.traj.t0

    return max(floor(Int, horizon / 0.4) + 2, 5)
end

using Clarabel # ECOS reports NUMERICAL_ERROR on a large fraction of these
               # problems, and sometimes returns an all NaN solution while
               # still reporting ALMOST_OPTIMAL. See docs/mpc-feasibility.md.
# using ECOS
# using Gurobi

import JuMP: MOI

"""
    lastUsableSolution(hist)

Rebuild an `SCPSolution` from the last PTR iteration whose conic subproblem
actually solved, or return `nothing` if there wasn't one.

`SCPToolbox` decides between `SCP_SOLVED` and `SCP_FAILED` by looking only at
the *last* subproblem, so one bad solve at the end discards every good iterate
before it — even though those are usually perfectly good trajectories. For a
controller that has to produce something every 0.25 s that is the wrong trade.

Also rejects non finite solutions: a subproblem can come back `ALMOST_OPTIMAL`
with a solution vector full of `NaN`, and `SCPToolbox` will happily
re-discretise about it (which is where the `SingularException`s come from).
"""
function lastUsableSolution(hist::SCPHistory)
    for subproblem in reverse(hist.subproblems)
        sol = subproblem.sol

        if (sol.status == MOI.OPTIMAL || sol.status == MOI.ALMOST_OPTIMAL) &&
           all(isfinite, sol.xd) && all(isfinite, sol.ud) && all(isfinite, sol.p)
            # SCPSolution takes the history and looks at its last entry, so hand
            # it a history truncated to this iteration.
            truncated = SCPHistory(hist.subproblems[1:subproblem.iter])

            return SCPSolution(truncated)
        end
    end

    return nothing
end

function ptr(mdl)
    # Problem definition
    pbm = TrajectoryProblem(mdl)
    define_problem!(pbm, :ptr)

    # Sized for the longest horizon the problem can pick, since the grid is
    # uniform in normalised time: if t_land stretches into the ballistic tail,
    # a grid sized for the burn alone would thin out over the powered phase too.
    N, Nsub = gridSize(mdl), 100 # dt can be set to 0.2 or even higher with little decrease in cost (velocity will only be a bit higher).
    # N can't be ≤ 1?
    iter_max = 50
    disc_method = FOH
    wvc, wtr = 5e3, 1e-2 # wtr is important, needs to be small but too small and we get problems.
    feas_tol = 1e-2
    ε_abs, ε_rel = 1e-5, 1e-3
    q_tr = Inf
    q_exit = Inf
    solver, options = Clarabel, Dict("verbose"=>false)
    # solver, options = ECOS, Dict("verbose"=>0)
    # MOI.Silent()
    # solver, options = Gurobi, Dict()#"OutputFlag"=>0)
    pars = PTR.Parameters(
        N, Nsub, iter_max, disc_method, wvc, wtr, ε_abs,
        ε_rel, feas_tol, q_tr, q_exit, solver, options)

    # Create and solve the problem
    ptr_pbm = PTR.create(pars, pbm)
    sol, hist = PTR.solve(ptr_pbm)

    if startswith(sol.status, string(SCP_FAILED))
        salvaged = lastUsableSolution(hist)

        if !isnothing(salvaged)
            return salvaged
        end
    end

    return sol
end

function scvx(mdl)
    # Problem definition
    pbm = TrajectoryProblem(mdl)
    define_problem!(pbm, :scvx)

    # PTR algorithm parameters
    N, Nsub = gridSize(mdl), 100
    iter_max = 50
    disc_method = FOH
    λ = 5e2
    ρ_0 = 0.0
    ρ_1 = 0.1
    ρ_2 = 0.7
    β_sh = 2.0
    β_gr = 2.0
    η_init = 1.0
    η_lb = 1e-8
    η_ub = 10.0
    feas_tol = 1e-2
    ε_abs, ε_rel = 1e-5, 1e-3
    q_tr = Inf
    q_exit = Inf
    solver, options = ECOS, Dict("verbose"=>0)
    pars = SCvx.Parameters(
        N, Nsub, iter_max, disc_method, λ, ρ_0, ρ_1, ρ_2, β_sh, β_gr, η_init,
        η_lb, η_ub, ε_abs, ε_rel, feas_tol, q_tr, q_exit, solver, options)

    # Create and solve the problem
    ptr_pbm = SCvx.create(pars, pbm)
    sol, history = SCvx.solve(ptr_pbm)
    return sol
end

function solveProblem(mdl = RocketProblem(), algo=:ptr)
    if algo == :ptr
        return ptr(mdl)
    elseif algo == :scvx
        return scvx(mdl)
    else
        error(algo, " not supported.")
    end
end