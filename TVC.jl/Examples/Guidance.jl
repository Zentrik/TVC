using TVC, SCPToolbox, Revise

# veh = RocketParameters(MomentArm = t -> [0; 0; -0.4]); # The parameters of the rocket
# atmos = Atmosphere(g = h -> [0; 0; -9.80655]);
veh = RocketParameters()
atmos = Atmosphere()
traj = RocketTrajectoryParameters()

mdl = RocketProblem(veh, atmos, traj)

sol = solveProblem(mdl);

printSolution(sol)
plotThrottle(sol)

sol.p

sample(sol.xc, 0 / veh.BurnTime)'

[traj.r0 + traj.v0 * 1. + atmos.g(traj.r0[3]) * 1^2/2; traj.v0 + atmos.g(traj.r0[3]) * 1]

traj = RocketTrajectoryParameters(r0 = [23.904477462043364, -6.928358096532524, 25.110234715614897], v0 = [3.8252602075575863, -2.8689451556681904, -9.7542178326703], q0 = [1.0, 0.0, 0.0, 0.0], ω0 = [0.0, 0.0, 0.0])
mdl = RocketProblem(veh, atmos, traj)

sol = solveProblem(mdl)

traj = RocketTrajectoryParameters(r0 = [23.904477462043364, -6.928358096532524, 25.110234715614897], v0 = [3.8252602075575863, -2.8689451556681904, -9.7542178326703], q0 = [0.9776157025914333, -0.1262391131789442, -0.1683188175719256, -5.423186447062239e-18], ω0 = [-0.446566778674036, -0.5954223715653811, 5.2138259556152435e-17])
mdl = RocketProblem(veh, atmos, traj)

sol = solveProblem(mdl)

printSolution(sol)
plotThrottle(sol)

using LinearAlgebra, SCPToolbox, Plots

t = LinRange(0, 1, 1000)
r = mapreduce(k -> sample(sol.xc, k)[veh.id_r], hcat, t)
Plots.plot(motorTime.(t, Ref(mdl)), r')

r = mapreduce(k -> sample(sol.uc, k)[veh.id_roll], hcat, t)
Plots.plot(motorTime.(t, Ref(mdl)), r')

sol.xd