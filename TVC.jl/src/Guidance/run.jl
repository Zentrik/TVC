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

export solveProblem

function ptr(mdl)
    # Problem definition
    pbm = TrajectoryProblem(mdl)
    define_problem!(pbm)

    N, Nsub = floor(Int, mdl.veh.BurnTime / 0.1) + 1, 15 # dt can be set to 0.2 or even higher with little decrease in cost (velocity will only be a bit higher).
    iter_max = 20
    disc_method = FOH
    wvc, wtr = 5e3, 1e-2 # wtr is important, needs to be small but too small and we get problems. 
    feas_tol = 10e-3
    ε_abs, ε_rel = 1e-5, 1e-3
    q_tr = Inf
    q_exit = Inf
    solver, options = ECOS, Dict("verbose"=>0)
    pars = PTR.Parameters(
        N, Nsub, iter_max, disc_method, wvc, wtr, ε_abs,
        ε_rel, feas_tol, q_tr, q_exit, solver, options)

    # Create and solve the problem
    ptr_pbm = PTR.create(pars, pbm)
    sol, history = PTR.solve(ptr_pbm)
    return sol
end

function solveProblem(mdl = RocketProblem())
    return ptr(mdl)
end