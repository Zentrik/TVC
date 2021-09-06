using Solvers
using Parser
using Utils

using LinearAlgebra
using ECOS
using PyPlot
using Colors
using Printf

include("./Rocket_Acceleration.jl")

export solve

function solve()
    # Define the problem
    pbm = TrajectoryProblem(nothing)
    problem_set_dims!(pbm, 6, 3, 2)

    g = [0; 0; -9.80655];

    # Initial guess
    r_0 = [20; -4; 30];	# position vector, m
    v_0 = [4; -3; 0];		# velocity vector, m/s
    r_N =[0; 0; 0];				# terminal position, m
    v_N =[0; 0; 0];				# terminal velocity, m
    problem_set_guess!(
        pbm, (N, pbm) -> begin
            # Cancel out gravity using, then solve for landing at t = 3.45
            t_coast = 0.0
            t_burn = 3.45
            tmp = [0; 0; - 2 * (r_0[3] + v_0[3] * t_burn) / t_burn^2] # [0.0, 0.0, -5.040957781978575]

            x = zeros(pbm.nx, N)
            u = zeros(pbm.nu, N)

            for k = 1:N
                t = (k - 1) / (N - 1) * t_burn
                x[:, k] = [r_0 + t * v_0 + t^2 * tmp / 2; v_0 + tmp * t]
                u[:, k] = -g + tmp # [0.0, 0.0, 4.765592218021425]
            end   

            p = [t_coast; t_burn]
            return x, u, p
        end)
    
    # Cost function
    # problem_set_running_cost!(
    #     pbm, :ptr, (t, k, x, u, p, pbm) -> #dot(x[4:6] - v_N, x[4:6] - v_N) # I don't think you can use norm, package needs to be updated or smth.
    #     #dot([zeros(3, 1); 2 * (x[4:6] - v_N)], f(x, u) * p[2]);
    #     #dot(x[4:6], g + u) * p[2]
    # )
    problem_set_terminal_cost!(
        pbm, (x, p, pbm) -> dot(x[4:6] - v_N, x[4:6] - v_N)
    )
    # Dynamics
    problem_set_dynamics!(
        pbm,
        # f
        (t, k, x, u, p, pbm) ->
            [x[4:6]; g + u]*p[2],
        # df/dx
        (t, k, x, u, p, pbm) ->
            [zeros(3, 3) I(3); zeros(3, 6)]*p[2],
        # df/du
        (t, k, x, u, p, pbm) ->
            [zeros(3, 3); I(3)]*p[2],
        # df/dp
        (t, k, x, u, p, pbm) ->
            [zeros(pbm.nx, 1) ones(pbm.nx, 1)])
    # Boundary conditions
    problem_set_bc!(
        pbm, :ic, # Initial condition
        (x, p, pbm) -> x - [r_0 + v_0 * p[1] + p[1]^2 * g / 2; v_0 + p[1] * g],
        (x, p, pbm) -> I(pbm.nx), # Jacobian wrt x
        (x, p, pbm) -> [[v_0 + g * p[1]; g] zeros(6, 1)] # Jacobian wrt p 
    )
    problem_set_bc!(
        pbm, :tc,  # Terminal condition
        (x, p, pbm) -> [x[3] - r_N[3]], # can't pass a scalar
        (x, p, pbm) -> [0 0 1 0 0 0]) 

    # Constraints s() <= 0
    # problem_set_s!(
    #     pbm, :ptr,
    #     # s
    #     (t, k, x, u, p, pbm) -> [-x[3]; -p[1]; p[1] - 5; 3.2 - p[2]; norm(u) - Acceleration(t * p[2])], # need to work out t coast max
    #     # ds/dx
    #     (t, k, x, u, p, pbm) -> [0 0 -1 0 0 0; zeros(4, pbm.nx)], 
    #     # ds/du
    #     (t, k, x, u, p, pbm) -> [zeros(4, pbm.nu); u' / norm(u)], # Symbolics.jacobian(norm(u), u)
    #     # ds/dp
    #     (t, k, x, u, p, pbm) -> [zeros(1, pbm.np); -1 0; 1 0; 0 -1; 0  t*central_fdm(5, 1)(Acceleration, t * p[2])] # df(t * p)/dp = df(z)/dz * dz/dp
    # )
    # problem_set_s!(
    #     pbm, :ptr,
    #     # s
    #     (t, k, x, u, p, pbm) -> [-x[3]; -p[1]; p[1] - 5; 3.2 - p[2]; norm(u) - 10], # need to work out t coast max
    #     # ds/dx
    #     (t, k, x, u, p, pbm) -> [0 0 -1 0 0 0; zeros(4, pbm.nx)], 
    #     # ds/du
    #     (t, k, x, u, p, pbm) -> [zeros(4, pbm.nu); u' / norm(u)], # Symbolics.jacobian(norm(u), u)
    #     # ds/dp
    #     (t, k, x, u, p, pbm) -> [zeros(1, pbm.np); -1 0; 1 0; 0 -1; 0 0] # df(t * p)/dp = df(z)/dz * dz/dp
    # )

    # Convex State Constraints
    problem_set_X!(
        pbm, (t, k, x, p, pbm, ocp) -> begin

        @add_constraint(
            ocp, NONPOS, "height >= 0", (x[3],), begin
            local height = arg[1]
            - height
            end)
        end)
    
    # Convex Input Constraints
    problem_set_U!(
        pbm, (t, k, u, p, pbm, ocp) -> begin

        @add_constraint(
            ocp, SOC, "Thrust Magnitude <= Max", (u,), begin # u is the descision variable for the constraint, but t isn't so don't include?
                local u = arg[1]
                [10.0; u]
            end)
        end)

    # problem_set_s!(
    #     pbm, :ptr,
    #     # s
    #     (t, k, x, u, p, pbm) -> [norm(u)^2 - Acceleration(t * p[2])^2], # need to work out t coast max
    #     # ds/dx
    #     (t, k, x, u, p, pbm) -> zeros(1, pbm.nx), 
    #     # ds/du
    #     (t, k, x, u, p, pbm) -> 2 * u',
    #     # ds/dp
    #     (t, k, x, u, p, pbm) -> [t * central_fdm(5, 1)(Acceleration, t * p[2])] # df(t * p)/dp = df(z)/dz * dz/dp
    # )    

    # Define the SCP algorithm parameters
    N, Nsub = floor(Int, 3.45 / 0.1) + 1, 10 # see LanderSolid.m, dt needs to be low. We need many degrees of freedom on u
    iter_max = 30
    disc_method = FOH
    wvc, wtr = 1e3, 0.11 # wtr is important, needs to be small but too small and we get problems
    feas_tol = 5e-3
    ε_abs, ε_rel = 1e-5, 1e-3
    q_tr = Inf
    q_exit = Inf
    solver, options = ECOS, Dict("verbose"=>0)
    pars = Solvers.PTR.Parameters(
        N, Nsub, iter_max, disc_method, wvc, wtr, ε_abs,
        ε_rel, feas_tol, q_tr, q_exit, solver, options)
    # Create and solve the problem
    ptr_pbm = Solvers.PTR.create(pars, pbm)
    sol, history = Solvers.PTR.solve(ptr_pbm)
    return sol
end