# TODO: Gimbal angular velocity constraint.
# Make quaternion discretisatin more accurate.

using Solvers
using Parser
using Utils

using LinearAlgebra
using ECOS

include("./Rocket_Acceleration.jl")
include("./Quaternions.jl")

# using Symbolics

using Plots
using NPZ

export solve, print, plot, save

r_0 = [20; -4; 30];	# position vector, m
v_0 = [4; -3; 0];		# velocity vector, m/s
q_0 = [1; 0; 0; 0];
w_0 = [0; 0; 0]

function solve(algo)
    pbm = TrajectoryProblem(nothing)
    problem_set_dims!(pbm, 13, 4, 1)
    set_scale!(pbm)
    set_integration!(pbm)

    g = [0; 0; -9.80655];

    r_N =[0; 0; 0];				# terminal position, m
    v_N =[0; 0; 0];				# terminal velocity, m
    q_N = [1; 0; 0; 0]; # or any roll of this.
    w_N = [0; 0; 0];

    problem_set_guess!(
        pbm, (N, pbm) -> begin
            t_coast = 0.0
            t_burn = 3.45

            # >> Parameter guess <<
            p = zeros(pbm.np)
            p[1] = t_coast

            # >> Input guess <<
            u = zeros(pbm.nu, N)
            for k = 1:N
                t = (k - 1) / (N - 1) * t_burn
                u[:, k] = [0; 0; Acceleration(t); 0]
            end

            # >> State guess <<
            x0 = [r_0 + v_0 * t_coast + g * t_coast^2/2; v_0 + g * t_coast; q_0; w_0];
            xf = [r_N; v_N; q_N; w_N]

            x = straightline_interpolate(x0, xf, N)

            # x = zeros(pbm.nx, N)
            # x[1:3, :] = straightline_interpolate(x0[1:3], xf[1:3], N)
            # x[4:6, :] = straightline_interpolate(x0[4:6], xf[4:6], N)
            
            # # for k = 1:N
            # #     mix = (k-1)/(N-1) # time normalised to 1
            # #     x[7:10, k] = slerp_quat(x0[7:10], xf[7:10], mix)
            # # end

            # x[7:10, :] = xf[7:10] .* ones(4, N)
            
            # x[11:13, :] = straightline_interpolate(x0[11:13], xf[11:13], N)

            return x, u, p 
        end
    )

    problem_set_terminal_cost!(
        pbm, (x, p, pbm) -> dot(x[4:6] - v_N, x[4:6] - v_N)
    )
    

    ## Testing Jacobians
    # using Symbolics
    # include("src/Quaternions.jl") 
    # Id = Diagonal([0.2, 0.2, 0.04])
    # invId = Diagonal([5.0, 5.0, 25.0])
    # g = [0; 0; -9.80655]

    # @variables r[1:3] v[1:3] u[1:4] quat[1:4] w[1:3];
    # r = Symbolics.scalarize(r);
    # v = Symbolics.scalarize(v);
    # quat = Symbolics.scalarize(quat);
    # w = Symbolics.scalarize(w);
    # u = Symbolics.scalarize(u);

    # x = [r; v; quat; w];
    
    # x_dot = [v; g + rotate(quat, u[1:3]); 1/2 * quatL(quat) * [0; w]; invId * (cross([0; 0; -0.5], u[1:3]) + [0; 0; u[4]] - cross(w, Id * w))]
    # A_sym = Symbolics.jacobian(x_dot, x)
    # B_sym = Symbolics.jacobian(x_dot, u)

    # f = build_function(x_dot, x, u, expression=Val{false})[1]
    # A = build_function(A_sym, x, u, expression=Val{false})[1]
    # B = build_function(B_sym, x, u, expression=Val{false})[1]   
    
    # function f_t(x, u)
    #     r = x[1:3]
    #     v = x[4:6]
    #     quat = x[7:10]
    #     w = x[11:13]
    #     return [v; g + rotate(quat, u[1:3]); 1/2 * quatL(quat) * [0; w]; inv(Id) * (cross([0; 0; -0.5], u[1:3]) + [0; 0; u[4]] - cross(w, Id * w))]
    # end

    # function A_t(x, u)
    #     r = x[1:3]
    #     v = x[4:6]
    #     quat = x[7:10]
    #     ω = x[11:13]

    #     A = zeros(13, 13)
    #     A[1:3, 4:6] = I(3)
    #     A[4:6, 7] = 2 * (quat[1] * u[1:3] + quat[2:4] × u[1:3])
    #     A[4:6, 8:10] = 2 * (quat[2:4]' * u[1:3] * I(3) + quat[2:4] * u[1:3]' - u[1:3] * quat[2:4]' - quat[1] * skew(u[1:3]))
    #     A[7:10, 7:10] = 1/2 * quatR([0; ω])
    #     A[7:10, 11:13] = 1/2 * quatL(quat)[:, 2:4]
    #     A[11:13, 11:13] = - indId * (skew(ω) * Id - skew(Id * ω))
    #     return A
    # end

    # function B_t(x, u)
    #     quat = x[7:10]

    #     B = zeros(13, 4)
    #     B[4:6, 1:3] = to_matrix(quat)
    #     B[11:13, 1:3] = inv(Id) * skew([0; 0; -0.5])
    #     B[11:13, 4] = inv(Id)[:, 3]
    #     return B
    # end

    # x = rand(13, 1);
    # u = rand(4, 1);
    # count(f(x, u) - f_t(x, u) .> 1e-10)
    # count(A(x, u) - A_t(x, u) .> 1e-10)
    # count(B(x, u) - B_t(x, u) .> 1e-10)

    Id = Diagonal([0.2, 0.2, 0.04])
    invId = Diagonal([5.0, 5.0, 25.0])
        
    # Dynamics
    problem_set_dynamics!(
        pbm,
        # f
        (t, k, x, u, p, pbm) -> begin
            r = x[1:3]
            v = x[4:6]
            quat = x[7:10]
            w = x[11:13]
            return [v; g + rotate(quat, u[1:3]); 1/2 * quatL(quat) * [0; w]; invId * (cross([0; 0; -0.5], u[1:3]) + [0; 0; u[4]] - cross(w, Id * w))] * 3.45
        end,
        # df/dx
        (t, k, x, u, p, pbm) -> begin
            r = x[1:3]
            v = x[4:6]
            quat = x[7:10]
            ω = x[11:13]

            A = zeros(pbm.nx, pbm.nx)
            A[1:3, 4:6] = I(3)
            A[4:6, 7] = 2 * (quat[1] * u[1:3] + quat[2:4] × u[1:3])
            A[4:6, 8:10] = 2 * (quat[2:4]' * u[1:3] * I(3) + quat[2:4] * u[1:3]' - u[1:3] * quat[2:4]' - quat[1] * skew(u[1:3]))
            A[7:10, 7:10] = 1/2 * quatR([0; ω])
            A[7:10, 11:13] = 1/2 * quatL(quat)[:, 2:4]
            A[11:13, 11:13] = - invId * (skew(ω) * Id - skew(Id * ω))
            A *= 3.45
            return A
        end,
        # df/du
        (t, k, x, u, p, pbm) -> begin
            quat = x[7:10]

            B = zeros(pbm.nx, pbm.nu)
            B[4:6, 1:3] = to_matrix(quat)
            B[11:13, 1:3] = invId * skew([0; 0; -0.5])
            B[11:13, 4] = invId[:, 3]
            B *= 3.45
            return B
        end,
        # df/dp
        (t, k, x, u, p, pbm) ->
            zeros(pbm.nx, pbm.np)*3.45)

    # Boundary conditions
    problem_set_bc!(
        pbm, :ic, # Initial condition
        (x, p, pbm) -> x - [r_0 + v_0 * p[1] + p[1]^2 * g / 2; v_0 + p[1] * g; q_0; w_0],
        (x, p, pbm) -> I(pbm.nx), # Jacobian wrt x
        (x, p, pbm) -> - [[v_0 + g * p[1]; g; zeros(pbm.nx - 6)] zeros(pbm.nx, pbm.np - 1)] # Jacobian wrt p 
    )

    problem_set_bc!(
        pbm, :tc,  # Terminal condition
        (x, p, pbm) -> [x[3] - r_N[3];
                        x[8:9] - q_N[2:3];
                        x[11:13] - w_N],
        (x, p, pbm) -> [zeros(1, 2) 1 zeros(1, pbm.nx - 3); 
                        zeros(2, 7) I(2) zeros(2, pbm.nx - 9);
                        zeros(3, 10) I(3)], 
        (x, p, pbm) -> zeros(6, pbm.np) # Jacobian wrt p
    ) 

    # Convex State Constraints
    problem_set_X!(
        pbm, (t, k, x, p, pbm, ocp) -> begin

        @add_constraint(
            ocp, NONPOS, "height >= 0", (x[3],), begin
            local height = arg[1]
            - height
            end)

        @add_constraint(
            ocp, NONPOS, "t_coast >= 0", (p[1],), begin
                local t_coast = arg[1]
                - t_coast
            end)
        
        @add_constraint( # shouldn't be necessary, but seems very very useful.
            ocp, SOC, "|quat| <= 1", (x[7:10],), begin
                local quat = arg[1]
                [1; quat]
            end)

        # @add_constraint(
        #     ocp, SOC, "|w| < w_max", (x[11:13],), begin
        #         local w = arg[1]
        #         [pi / 2; w]
        #     end)
        end)

    # Convex Input Constraints
    problem_set_U!(
        pbm, (t, k, u, p, pbm, ocp) -> begin

        @add_constraint(
            ocp, SOC, "Thrust Magnitude <= Max", (u[1:3],), begin # u is the descision variable for the constraint, but t isn't so don't include?
                local u = arg[1]
                # [12; u]
                [Acceleration(t * 3.45); u]
            end)

        @add_constraint(
            ocp, L1, "|roll torque| < max", (u[4],), begin
                local u = arg[1]
                [0.1; u]
            end)
            
        @add_constraint(
            ocp, SOC, "Thrust Gimal angle < delta_max", (u[1:3],), begin
                local u = arg[1]
                [u[3] / cos(pi * 5/180); u]
            end)
        end) 

    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # :: SCvx algorithm parameters :::::::::::::::::::::::::::::::::::::::::::::::::
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    if algo == :scvx
        N = floor(Int, 3.45 / 0.1) + 1
        Nsub = 100
        iter_max = 30
        disc_method = FOH
        λ = 5e3
        ρ_0 = 0.0
        ρ_1 = 0.1
        ρ_2 = 0.7
        β_sh = 2.0
        β_gr = 2.0
        η_init = 1.0
        η_lb = 1e-8
        η_ub = 10.0
        ε_abs = 1e-8
        ε_rel = 1e-5
        feas_tol = 5e-3
        q_tr = Inf
        q_exit = Inf
        solver = ECOS
        solver_options = Dict("verbose"=>0, "maxit"=>1000)
        pars = Solvers.SCvx.Parameters(N, Nsub, iter_max, disc_method, λ, ρ_0, ρ_1, ρ_2, β_sh, β_gr,
                            η_init, η_lb, η_ub, ε_abs, ε_rel, feas_tol, q_tr,
                            q_exit, solver, solver_options)

        # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        # :: Solve trajectory generation problem ::::::::::::::::::::::::::::::::::::::
        # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        scvx_pbm = Solvers.SCvx.create(pars, pbm)
        sol, history = Solvers.SCvx.solve(scvx_pbm)
    elseif algo == :ptr
        N, Nsub = floor(Int, 3.45 / 0.4) + 1, 15 # see LanderSolid.m, dt needs to be low. We need many degrees of freedom on u. We set dt high so tvc_dot is low. # Nsub is how many points to use to calculate discretization between each timestep
        iter_max = 20
        disc_method = FOH
        wvc, wtr = 5e3, 1e-2 # wtr is important, needs to be small but too small and we get problems. 
        feas_tol = 10e-3
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
    end
    return sol
end

function set_scale!(pbm::TrajectoryProblem)::Nothing #VERY IMPORTANT
    advise! = problem_advise_scale!

    # States
    advise!(pbm, :state, 1, (-1000.0, 1000.0))
    advise!(pbm, :state, 2, (-1000.0, 1000.0))
    advise!(pbm, :state, 3, (0.0, 1000.0))
    advise!(pbm, :state, 4, (-100.0, 100.0))
    advise!(pbm, :state, 5, (-100.0, 100.0))
    advise!(pbm, :state, 6, (-100.0, 100.0))

    advise!(pbm, :state, 7, (-1.0, 1.0))
    advise!(pbm, :state, 8, (-1.0, 1.0))
    advise!(pbm, :state, 9, (-1.0, 1.0))
    advise!(pbm, :state, 10, (-1.0, 1.0))
    advise!(pbm, :state, 11, (-10.0, 10.0))
    advise!(pbm, :state, 12, (-10.0, 10.0))
    advise!(pbm, :state, 13, (-10.0, 00.0))

    # Inputs
    advise!(pbm, :input, 1, (-25.0, 25.0))
    advise!(pbm, :input, 2, (-25.0, 25.0))
    advise!(pbm, :input, 3, (-25.0, 25.0))
    advise!(pbm, :input, 4, (-5.0, 5.0))

    # Parameters
    advise!(pbm, :parameter, 1, (0.0, 10.0))

    return nothing
end

function set_integration!(pbm::TrajectoryProblem)::Nothing

    # Quaternion re-normalization on numerical integration step
    problem_set_integration_action!( # i think this only affects .xc
        pbm, 7:10,
        (q, pbm) -> begin
            qn = q/norm(q)
            return qn
        end)

    return nothing
end # function


function print(solution)
    println("Coast time (s): ", solution.p[1])

    v_N =[0; 0; 0];
    println("Impact Velocity Magnitude (m/s): ", solution.cost^0.5)
end

function plot(solution)
    t = solution.td * 3.45
    Plots.plot(t, ThrottleLevel(norm.(eachcol(solution.ud[1:3, :])), Acceleration(t) ), title = "Throttle Level" )
end

function save(solution)
    npzwrite("x.npy", transpose([[r_0; v_0; q_0; w_0] solution.xd]))
    npzwrite("u.npy", transpose([[0; 0; 0; 0] solution.ud]))
    npzwrite("t.npy", transpose([0; solution.p[1] .+ 3.45 * solution.td]))
end