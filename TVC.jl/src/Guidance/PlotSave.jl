using SCPToolbox
using LinearAlgebra

using Plots
using NPZ
using DelimitedFiles

import TVC: motorTime

export printSolution, plotThrottle, saveTrajectory

"""
    motorTime(t, mdl)

Gets time since rocket motor was ignited.
"""

function motorTime(t, mdl::RocketProblem)
    return mdl.traj.t0 + t * (mdl.veh.BurnTime - mdl.traj.t0)
end

function printSolution(solution)
    println("Coast time (s): ", solution.p[1])

    println("Impact Velocity Magnitude (m/s): ", solution.cost^0.5)

    println("True Impact Velocity Magnitude (m/s): ", norm(sample(solution.xc, solution.xc.t[end])[4:6]) ) # solution.xc.t[end] is just 1.0, solution.xc[3] at end isn't 0.
end

function plotThrottle(solution, mdl::RocketProblem = RocketProblem())
    Plots.plot(motorTime.(solution.td, Ref(mdl)), norm.(eachcol(solution.xd[mdl.veh.id_T, :])), title = "Throttle Level" )

    t = motorTime.(LinRange(0, 1, 1000), Ref(mdl))
    Plots.plot!(t, [norm(sample(solution.xc, k)[mdl.veh.id_T]) for k in LinRange(0, 1, 1000)], title="Throttle Level")

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

function saveTrajectory(solution, mdl::RocketProblem = RocketProblem())
    veh = mdl.veh
    traj = mdl.traj

    N = size(solution.xd)[2]

    x0 = zeros(13)
    x0[veh.id_r] .= traj.r0
    x0[veh.id_v] .= traj.v0
    x0[veh.id_quat] .= traj.q0
    x0[veh.id_ω] = traj.ω0


    npzwrite("x.npy", transpose([x0 solution.xd[vcat(veh.id_r, veh.id_v, veh.id_quat, veh.id_ω), :]]) .- [solution.xd[veh.id_r, end]' .* ones(N + 1, 3) zeros(N + 1, 10)]) # set final position to 0, so rocket lands on pad.
    npzwrite("u.npy", [T_0' * 0; solution.xd[veh.id_T, :]' .* veh.Thrust(motorTime.(solution.td, Ref(mdl)))] )
    npzwrite("t.npy", transpose([0; solution.p[veh.id_tcoast] .+ motorTime.(solution.td, Ref(mdl))]))

    # Save for simulink 

    t_fine = 0:0.01:veh.BurnTime
    x_fine = transpose(reduce(hcat, sample.(Ref(solution.xc), t_fine / max(t_fine...)))) # so we can do linear interpolation in matlab and be accurate
    u_fine = transpose(reduce(hcat, sample.(Ref(solution.uc), t_fine / max(t_fine...))))
    writedlm("x.csv", hcat(t_fine, x_fine[:, vcat(veh.id_r, veh.id_v, veh.id_quat, veh.id_ω)]), ',')

    Moments = zeros(length(t_fine), 3)
    for (k, t) in enumerate(t_fine)
        Moments[k, :] = cross(veh.MomentArm(t), x_fine[k, veh.id_T]) * veh.Thrust(t) + [0; 0; u_fine[k, veh.id_roll]]
    end
    writedlm("moments.csv", hcat(t_fine, Moments), ',')

    writedlm("u.csv", hcat([0; solution.p[veh.id_tcoast] - 1e-10; solution.p[veh.id_tcoast] .+ motorTime.(solution.td, Ref(mdl))], [[0 0 0 0]; [0 0 0 0]; solution.ud']), ',')
end