using Solvers
using Parser
using Utils

using LinearAlgebra
using ECOS
using PyPlot
using Colors
using Printf

export solve

_ro = 0.35
_co = [-0.1; 1]
_carl, _carw = 0.2, 0.1

function solve()
    # Define the problem
    pbm = TrajectoryProblem(nothing)
    problem_set_dims!(pbm, 3, 2, 1)
    # Initial guess
    _x0 = [0; 0; 0] # Initial state
    _xf = [0; 2; 0] # Terminal state
    problem_set_guess!(
        pbm, (N, pbm) -> begin
            pars = pbm.mdl
            x = straightline_interpolate(_x0, _xf, N)
            idle = zeros(pbm.nu)
            u = straightline_interpolate(idle, idle, N)
            p = zeros(pbm.np)
            return x, u, p
        end)
    # Cost function
    problem_set_running_cost!(
        pbm, :ptr, (t, k, x, u, p, pbm) -> u'*u)
    # Dynamics
    _tf = 3
    problem_set_dynamics!(
        pbm,
        # f
        (t, k, x, u, p, pbm) ->
            [u[1]*sin(x[3]);
             u[1]*cos(x[3]);
             u[2]]*_tf,
        # df/dx
        (t, k, x, u, p, pbm) ->
            [0 0 u[1]*cos(x[3]);
             0 0 -u[1]*sin(x[3]);
             0 0 0]*_tf,
        # df/du
        (t, k, x, u, p, pbm) ->
            [sin(x[3]) 0;
             cos(x[3]) 0;
             0 1]*_tf,
        # df/dp
        (t, k, x, u, p, pbm) ->
            zeros(pbm.nx, pbm.np))
    # Obstacle constraint
    problem_set_s!(
        pbm, :ptr,
        # s
        (t, k, x, u, p, pbm) -> [(_ro+_carw/2)^2-(x[1]-_co[1])^2-(x[2]-_co[2])^2],
        # ds/dx
        (t, k, x, u, p, pbm) -> collect([-2*(x[1]-_co[1]); -2*(x[2]-_co[2]); 0]'),
        # ds/du
        (t, k, x, u, p, pbm) -> zeros(1, pbm.nu),
        # ds/dp
        (t, k, x, u, p, pbm) -> zeros(1, pbm.np))
    # Boundary conditions
    problem_set_bc!(
        pbm, :ic, # Initial condition
        (x, p, pbm) -> x-_x0,
        (x, p, pbm) -> I(pbm.nx))
    problem_set_bc!(
        pbm, :tc,  # Terminal condition
        (x, p, pbm) -> x-_xf,
        (x, p, pbm) -> I(pbm.nx))
    # Define the SCP algorithm parameters
    N, Nsub = 11, 10
    iter_max = 30
    disc_method = FOH
    wvc, wtr = 1e3, 1e0
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

function plot(sol)
    # Trajectory plot
    ctres, overlap = 1000, 3
    N = size(sol.xd, 2)
    xct = hcat([sample(sol.xc, t) for t in LinRange(0, 1, ctres)]...)
    vct = vcat([sample(sol.uc, t)[1] for t in LinRange(0, 1, ctres)]...)
    cmap = generate_colormap("inferno"; minval=minimum(vct), maxval=maximum(vct))
    plot_options = Dict("xtick.labelsize"=>9,
                        "ytick.labelsize"=>9,
                        "axes.labelsize"=>11)
    fig = create_figure((4, 4); options = plot_options)
    ax = setup_axis!(111, xlabel="\$x\$ [m]", ylabel="\$y\$ [m]",
                     axis="equal", cbar=cmap, clabel="Velocity, \$v\$ [m/s]",
                     cbar_aspect=40)
    ax.plot(sol.xd[1, :], sol.xd[2, :],
            linestyle="none", marker="o", markerfacecolor=DarkBlue,
            markeredgecolor="white", markeredgewidth=0.2, markersize=3,
            zorder=20)
    line_segs = Vector{Matrix}(undef, 0)
    line_clrs = Vector{NTuple{4, Real}}(undef, 0)
    for k=1:ctres-overlap
        push!(line_segs, xct[1:2, k:k+overlap]')
        push!(line_clrs, cmap.to_rgba(vct[k]))
    end
    trajectory = PyPlot.matplotlib.collections.LineCollection(
        line_segs, zorder=10, colors = line_clrs, linewidths=3,
        capstyle="round")
    ax.add_collection(trajectory)
    Rect = PyPlot.matplotlib.patches.Rectangle
    for k=1:N
        local xl, xw = [1;1;-1;-1;1]*_carl/2, [1;-1;-1;1;1]*_carw/2
        local yl, yw = [1;1;-1;-1;1]*_carl/2, [-1;1;1;-1;-1]*_carw/2
        local ang = sol.xd[3,k]
        local xc = sol.xd[1,k].+xl.*sin(ang).+xw.*cos(ang)
        local yc = sol.xd[2,k].+yl.*cos(ang).+yw.*sin(ang)
        ax.fill(xc, yc,
                linewidth=1,
                edgecolor=DarkBlue,
                facecolor=rgb2pyplot(parse(RGB, Red), a=0.5),
                zorder=6)
    end
    ang = LinRange(0, 2*pi, 100)
    obs = ([cos.(ang)'; sin.(ang)']*_ro).+_co
    ax.fill(obs[1, :], obs[2, :],
            linewidth=1,
            edgecolor=Blue,
            facecolor=rgb2pyplot(parse(RGB, Green), a=0.5),
            zorder=5)
    save_figure("../figures/dubin_trajectory.pdf")

    return nothing
end
