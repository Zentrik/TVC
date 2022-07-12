using TVC, SCPToolbox

# veh = RocketParameters(MomentArm = t -> [0; 0; -0.4]); # The parameters of the rocket
# atmos = Atmosphere(g = h -> [0; 0; -9.80655]);
veh = RocketParameters()
atmos = Atmosphere()
traj = RocketTrajectoryParameters()

mdl = RocketProblem(veh, atmos, traj)

sol = solveProblem(mdl);

printSolution(sol)
plotThrottle(sol)

using LinearAlgebra, SCPToolbox, Plots

t = LinRange(0, 1, 1000)
r = mapreduce(k -> sample(sol.xc, k)[veh.id_r], hcat, t)
Plots.plot(motorTime.(t, Ref(mdl)), r')

r = mapreduce(k -> sample(sol.uc, k)[veh.id_roll], hcat, t)
Plots.plot(motorTime.(t, Ref(mdl)), r')

sol.xd