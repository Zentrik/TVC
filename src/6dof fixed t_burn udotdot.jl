# TODO: Make quaternion discretisation more accurate. Without constraint on quaternion norm, the norm of the quaternion strays from 1, I think this is due to the Jacobians of dynamics being evaluated about a reference trajectory that has non unit quaternions.

using Solvers
using Parser
using Utils

using LinearAlgebra
using ECOS

include("./Rocket_Acceleration.jl")
include("./Quaternions.jl")

using Plots
using NPZ

export solve, print, plot, save

r_0 = [20; -4; 30];	# position vector, m
v_0 = [4; -3; 0];		# velocity vector, m/s
q_0 = [1; 0; 0; 0];
w_0 = [0; 0; 0]
T_0 = [0; 0; 1] # acceleration from Thrust normalised
Tdot_0 = [0; 0; 0]

# Per Successive Convexification for Mars 6-DoF Powered Descent Landing Guidance, 2017. Set control to second derivative of thrust vector
# This allows me to add a constraint on the tvc gimbal rate and it adds more degrees of freedom on the control/ allows for more complex control for a given N (no. of discretisation steps) improving cost.


function solve(algo)
    pbm = TrajectoryProblem(nothing)
    problem_set_dims!(pbm, 19, 4, 1)
    set_scale!(pbm)

    g = [0; 0; -9.80655];

    r_N =[0; 0; 0];				# terminal position, m
    v_N =[0; 0; 0];				# terminal velocity, m
    q_N = [1; 0; 0; 0]; # or any roll of this.
    w_N = [0; 0; 0];

    T_N = [0; 0; 1]
    Tdot_0 = [0; 0; 0]

    problem_set_guess!(
        pbm, (N, pbm) -> begin
            global dt = 3.45 / (N - 1) # for convex constraints

            t_coast = 0.0
            t_burn = 3.45

            x0 = [r_0 + v_0 * t_coast + g * t_coast^2/2; v_0 + g * t_coast; q_0; w_0; T_0; Tdot_0];
            xf = [r_N; v_N; q_N; w_N; T_N; Tdot_0]
            
            x = zeros(pbm.nx, N)
            u = zeros(pbm.nu, N)
            
            p = [t_coast]

            for k = 1:N
                t = (k - 1) / (N - 1) * t_burn
                x[:, k] = (t_burn - t) / t_burn * x0 + t / t_burn * xf
                u[:, k] = [0; 0; 0; 0]
            end

            return x, u, p 
        end
    )

    problem_set_terminal_cost!(
        pbm, (x, p, pbm) -> dot(x[4:6] - v_N, x[4:6] - v_N)
    )
    
    Id = Diagonal([0.2, 0.2, 0.04])
    invId = Diagonal([5.0, 5.0, 25.0])
    g = [0; 0; -9.80655]

    ## Testing Jacobians
    # using Symbolics, LinearAlgebra
    # include("src/Quaternions.jl") 

    # @variables r[1:3] v[1:3] u[1:4] quat[1:4] w[1:3] T[1:3] T_dot[1:3] Thrust;
    # r = Symbolics.scalarize(r);
    # v = Symbolics.scalarize(v);
    # quat = Symbolics.scalarize(quat);
    # w = Symbolics.scalarize(w);
    # T = Symbolics.scalarize(T);
    # T_dot = Symbolics.scalarize(T_dot);

    # u = Symbolics.scalarize(u);

    # x = [r; v; quat; w; T; T_dot];
    
    # x_dot = [v; g + rotate(quat, T) * Thrust; 1/2 * quatL(quat) * [0; w]; invId * (cross([0; 0; -0.5], T * Thrust) + [0; 0; u[4]] - cross(w, Id * w)); T_dot; u[1:3]]
    # A_sym = Symbolics.jacobian(x_dot, x)
    # B_sym = Symbolics.jacobian(x_dot, u)

    # f = build_function(x_dot, x, u, Thrust, expression=Val{false})[1]
    # A = build_function(A_sym, x, u, Thrust, expression=Val{false})[1]
    # B = build_function(B_sym, x, u, Thrust, expression=Val{false})[1]   
    
    function f_t(x, u, Thrust)
        r = x[1:3]
        v = x[4:6]
        quat = x[7:10]
        w = x[11:13]
        T = x[14:16]
        T_dot = x[17:19]

        return [v; g + rotate(quat, T) * Thrust; 1/2 * quatL(quat) * [0; w]; invId * (cross([0; 0; -0.5], T * Thrust) + [0; 0; u[4]] - cross(w, Id * w)); T_dot; u[1:3]]
    end

    function A_t(x, u, Thrust)
        r = x[1:3]
        v = x[4:6]
        quat = x[7:10]
        ω = x[11:13]
        T = x[14:16]
        T_dot = x[17:19]

        A = zeros(19, 19)
        A[1:3, 4:6] = I(3)
        A[4:6, 7] = Thrust * 2 * (quat[1] * T + quat[2:4] × T)
        A[4:6, 8:10] = Thrust * 2 * (quat[2:4]' * T * I(3) + quat[2:4] * T' - T * quat[2:4]' - quat[1] * skew(T))
        A[4:6, 14:16] = Thrust * to_matrix(quat)
        A[7:10, 7:10] = 1/2 * quatR([0; ω])
        A[7:10, 11:13] = 1/2 * quatL(quat)[:, 2:4]
        A[11:13, 11:13] = - invId * (skew(ω) * Id - skew(Id * ω))
        A[11:13, 14:16] = invId * skew([0; 0; -0.5]) * Thrust
        A[14:16, 17:19] = I(3)

        return A
    end

    function B_t(x, u, Thrust)
        B = zeros(19, 4)
        B[11:13, 4] = invId[:, 3]
        B[17:19, 1:3] = I(3)

        return B
    end

    # x = rand(19, 1);
    # u = rand(4, 1);
    # Thrust = rand()
    # count(f(x, u, Thrust) - f_t(x, u, Thrust) .> 1e-10)
    # count(A(x, u, Thrust) - A_t(x, u, Thrust) .> 1e-10)
    # count(B(x, u, Thrust) - B_t(x, u, Thrust) .> 1e-10)

    # Dynamics
    problem_set_dynamics!(
        pbm,
        # f
        (t, k, x, u, p, pbm) -> f_t(x, u, Acceleration(t * 3.45)) * 3.45,
        # df/dx
        (t, k, x, u, p, pbm) -> A_t(x, u, Acceleration(t * 3.45)) * 3.45,
        # df/du
        (t, k, x, u, p, pbm) -> B_t(x, u, Acceleration(t * 3.45)) * 3.45,
        # df/dp
        (t, k, x, u, p, pbm) ->
            zeros(pbm.nx, pbm.np)*3.45)

     # Quaternion re-normalization on numerical integration step
     problem_set_integration_action!(
        pbm, 7:10,
        (q, pbm) -> begin
            qn = q/norm(q)
            return qn
        end)

    # Boundary conditions
    problem_set_bc!(
        pbm, :ic, # Initial condition
        (x, p, pbm) -> x - [r_0 + v_0 * p[1] + p[1]^2 * g / 2; v_0 + p[1] * g; q_0; w_0; T_0; Tdot_0],
        (x, p, pbm) -> I(pbm.nx), # Jacobian wrt x
        (x, p, pbm) -> - [[v_0 + g * p[1]; g; zeros(pbm.nx - 6)] zeros(pbm.nx, pbm.np - 1)] # Jacobian wrt p 
    )

    problem_set_bc!(
        pbm, :tc,  # Terminal condition
        (x, p, pbm) -> [x[3] - r_N[3];
                        x[8:9] - q_N[2:3];
                        x[11:13] - w_N;
                        x[14:16] - T_N],
        (x, p, pbm) -> begin
            J = zeros(9, pbm.nx)
            J[1, 3] = 1
            J[2:3, 8:9] = I(2)
            J[4:6, 11:13] = I(3)
            J[7:9, 14:16] = I(3)
            return J
        end,
        (x, p, pbm) -> zeros(9, pbm.np) # Jacobian wrt p
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

        @add_constraint(
            ocp, SOC, "Thrust Magnitude <= Max", (x[14:16],), begin # we have say thrust <= 1, as we want it normalised
                local Thrust = arg[1]
                [1; Thrust]
            end)

        @add_constraint(
            ocp, SOC, "Thrust Gimal angle < delta_max", (x[14:16],), begin
                local Thrust = arg[1]
                [Thrust[3] / cos(pi * 5/180); Thrust]
            end)
        
        @add_constraint( # shouldn't be necessary, but seems very very useful.
            ocp, SOC, "|quat| <= 1", (x[7:10],), begin
                local quat = arg[1]
                [1; quat]
            end)

        @add_constraint(
        ocp, SOC, "TVC angular velocity <= Max", (x[17:19],), begin
            local u = arg[1]
            [deg2rad(5); u] # Angular velocity is r × v / ||r||², v = u, assume u is ⊥ r and ||r||² = 1, so angular velocity is v = u. 
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
            ocp, L1, "|roll torque| < max", (u[4],), begin
                local u = arg[1]
                [0.1; u]
            end)

        @add_constraint(
            ocp, SOC, "TVC Acceleration <= Max", (u[1:3],), begin
                local u = arg[1]
                [deg2rad(10); u] # Angular Acceleration is r × a / ||r||², a = u, assume u is ⊥ r and ||r||² = 1, so angular acceleration is a = u. 
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
        N, Nsub = floor(Int, 3.45 / 0.1) + 1, 15 # see LanderSolid.m, dt needs to be low. We need many degrees of freedom on u. We set dt high so tvc_dot is low. # Nsub is how many points to use to calculate discretization between each timestep
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

    return nothing
end

function print(solution)
    println("Coast time (s): ", solution.p[1])

    v_N = [0; 0; 0];
    println("Impact Velocity Magnitude (m/s): ", solution.cost^0.5)
end

function plot(solution)
    t_burn = 3.45
    t = solution.td * t_burn
    Plots.plot(t, norm.(eachcol(solution.xd[14:16, :])), title = "Throttle Level" )

    t = LinRange(0, 1, 1000) * t_burn
    Plots.plot!(t, [norm(sample(solution.xc, k)[14:16]) for k in t / t_burn], title="Throttle Level")

    # Plots.plot(t, norm.(eachcol(solution.xd[17:19, :])) ./ dt)

    # ang_vel = [norm(cross(solution.xd[14:16, i], solution.xd[17:19, i]) / norm(solution.xd[14:16, i])^2) for i = 1:size(solution.ud)[2]]
    # Plots.plot(solution.td .* t_burn, rad2deg.(ang_vel), title = "TVC Velocity (degrees)")


    # Gives same result as above after bugfixes

    # N = length(solution.td)
    # tvc_dot_max = similar(solution.td[1:N-1])
    # dt = solution.td[2] * 3.45

    # for k in eachindex(tvc_dot_max[1:N - 1])
    #     tvc_dot_max[k] = atand(norm(solution.xd[14:16, k] × solution.xd[14:16, k + 1]), solution.xd[14:16, k] ⋅ solution.xd[14:16, k + 1]) / dt # https://discourse.julialang.org/t/how-to-find-angle-between-two-vectors/68149/3
    #     # acosd(u_normalised[1:3, k] ⋅ u_normalised[1:3, k + 1] ) / dt
    # end

    # Plots.plot(solution.td[1:N - 1] * t_burn, tvc_dot_max)

    # Plots.plot(t, [rad2deg(norm(sample(solution.uc, k)[1:3])) for k in t / t_burn], title = "TVC Acceleration")
end

function save(solution)
    N = size(solution.xd)[2]
    npzwrite("x.npy", transpose([[r_0; v_0; q_0; w_0] solution.xd[1:13, :]]) .- [solution.xd[1:3, end]' .* ones(N + 1, 3) zeros(N + 1, 10)]) # set final position to 0, so rocket lands on pad.
    npzwrite("u.npy", [T_0' * 0; solution.xd[14:16, :]' .* Acceleration(solution.td * 3.45)] )
    npzwrite("t.npy", transpose([0; solution.p[1] .+ 3.45 * solution.td]))
end