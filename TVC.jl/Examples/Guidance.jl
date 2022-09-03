using TVC, SCPToolbox

veh = RocketParameters()
atmos = Atmosphere()
# atmos = Atmosphere(g = h -> [0; 0; -9.80655]); # Let g be a constant.
traj = RocketTrajectoryParameters()
# traj = RocketTrajectoryParameters(r0=[29.748114128035194, -11.314733910010972, 8.35342199993651], v0=[3.0171645950297212, -2.2567065996045876, -7.7594061952317785], q0=[0.9944172311103161, -0.06223739713880662, -0.07377484414430403, 0.042639761140275734], ω0=[-0.048223408675021756, -0.0545934384189852, 0.0365528711510246], T0=[-0.008838309484347663, 0.007928331573572896, 0.9999292376990102], Ṫ0=[0.0024214667602295878, -0.005573117719239776, 6.566557159870377e-5], t0=1.174078871814609, MotorFired=true);
# traj = RocketTrajectoryParameters(r0=[30.432078228827013, -11.824638836538547, 6.532740431280124], v0=[2.4624700463201536, -1.827961722129504, -6.816902726306204], q0=[0.9955741335293473, -0.05767366213382853, -0.07405014956352449, 0.004740114272751105], ω0=[0.06144204383116665, 0.07297527336306442, -0.35178497728103136], T0=[-0.0037608375882585607, 0.0037602884190028668, 0.9999231804329158], Ṫ0=[0.0333300182491607, -0.02758081010211138, 0.0015782812806192557], t0=1.424078871814609, MotorFired=true);
# traj = RocketTrajectoryParameters(r0=[31.559679415909518, -12.744739025910384, 3.630009892906172], v0=[1.658040221310699, -1.2618396720481064, -4.88140052568886], q0=[0.952461545091973, -0.0635383078438689, -0.04037860870066029, 0.2952108678838691], ω0=[0.08642184793974272, -0.005138395683974999, 0.12974083410017367], T0=[0.0040260160903892684, -0.003209297213079192, 0.9999492554481306], Ṫ0=[-0.005473940875827542, -0.01754906448034367, 0.0004519002932869254], t0=1.9236975910646925, MotorFired=true);
traj = RocketTrajectoryParameters(r0=[31.566947396947118, -12.713656291133669, 3.6300670318646358], v0=[1.6394448466332443, -1.2400501440115568, -4.883350399254044], q0=[0.9703116421702913, -0.05961284105483699, -0.047929027450857774, 0.2294437503567599], ω0=[0.08301297381039931, 0.014859974244157158, 0.13787449867442705], T0=[0.0029974169890587643, -0.002516070632471221, 0.9999977467479241], Ṫ0=[-0.0007944251430647212, -0.01197456611320066, -4.1720144142818e-5], t0=1.923747607139962, MotorFired=true);
traj = RocketTrajectoryParameters(r0=[31.56853948104539, -12.67731326024456, 3.628176975346708], v0=[1.6315425568554713, -1.218777300239065, -4.886976236334344], q0=[0.9973203826563598, -0.04377278034541103, -0.058617375031693286, 1.9463043707261873e-5], ω0=[0.060611729424207804, 0.07500123347971271, 0.0006343190100529849], T0=[0.0039035537230018907, -0.0030444017984405662, 0.9999878268814247], Ṫ0=[0.010911557712167579, -0.009233178037154224, -5.88577184501331e-5], t0=1.923747607139962, MotorFired=true);

mdl = RocketProblem(veh, atmos, traj)

solguidance = solveProblem(mdl);
printSolution(solguidance)
plotThrottle(solguidance)

using Plots

# t = LinRange(0, 1, 1000)
# r = mapreduce(k -> sample(sol.xc, k)[veh.id_r], hcat, t)
# Plots.plot(motorTime.(t, Ref(mdl)), r')

# r = mapreduce(k -> sample(sol.uc, k)[veh.id_roll], hcat, t)
# Plots.plot(motorTime.(t, Ref(mdl)), r')

Plots.plot!(traj.t0 .+ (mdl.veh.BurnTime - traj.t0) * solguidance.td, solguidance.xd[mdl.veh.id_r, :]')

# Plots.plot!(sol.prob.p.MotorIgnitionTime .+ traj.t0 .+ (mdl.veh.BurnTime - traj.t0) * solguidance.td, solguidance.xd[mdl.veh.id_v, :]')