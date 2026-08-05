using TVC, DifferentialEquations, Plots, LinearAlgebra

using Parameters

@with_kw struct ODEParametersFeas{R, V}
    veh::RocketParameters = RocketParameters()
    atmos::Atmosphere = Atmosphere()
    traj::RocketTrajectoryParameters = RocketTrajectoryParameters();

    Aero::Bool = true
    wind::V = zeros(3)

    ground::Bool = false

    MotorIgnitionTime::R = 0.0
    Control = (x, p, t) -> (force=[0; 0; veh.Thrust(t)], torque=zeros(3))
end

#   Can we land on the ground at end of burn.
#   ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡

traj = RocketTrajectoryParameters(ω0=[0, 0, 0.]);

traj = RocketTrajectoryParameters(r0=[22.753125968073366, -6.064844476055015, 27.600343655405403], v0=[3.859682362140023, -2.8947617716050167, -6.8455376409576365], q0=[0.9942734108276859, -0.0641197194907359, -0.08549295932098121, 1.481998974917953e-18], ω0=[-0.3699397941024044, -0.49325305880320575, -2.5565605675535743e-17], t0=0.0, MotorFired=false);

traj=RocketTrajectoryParameters(r0=[22.366130464744867, -5.774597848558639, 28.23627501803703], v0=[3.8806161631580847, -2.910462122368562, -5.872678737716989], q0=[0.9969377999299179, -0.046919167783644675, -0.06255889037819291, -7.544599728408613e-20], ω0=[-0.3199687185611176, -0.4266249580814899, 2.3758547531861336e-17])

traj=RocketTrajectoryParameters(r0=[24.989955030524015, -7.742466272892784, 22.023941100697545], v0=[3.880616164336729, -2.9104621217169893, -12.503308539930842], q0=[0.9969377604325815, -0.04691916581432842, -0.06255888701646095, -2.491975514828937e-10], ω0=[-0.3199687174763781, -0.42662495386317545, 2.504663143554353e-13])

traj=RocketTrajectoryParameters(r0=[26.782390393208537, -8.851442025590758, 8.149279787392047], v0=[-0.1435211379066393, 0.5706078127399808, -10.186068497061449], q0=[-0.08249345902837653, -0.024647815028127296, 0.02246904983541904, -0.9960333610031583], ω0=[-0.5726778952764043, -0.8717586111590575, -5.864091233219142], t0=1.3495278309691208, MotorFired=true);

x₀ = [traj.r0; traj.v0; traj.q0; traj.ω0][[veh.id_r; veh.id_v; veh.id_quat; veh.id_ω]];

p = ODEParametersFeas(traj=traj)

#   Solve ODE
#   ≡≡≡≡≡≡≡≡≡≡≡

tspan = (traj.t0, veh.BurnTime)

prob = ODEProblem(f!, x₀, tspan, p)

sol = DifferentialEquations.solve(prob, reltol=1e-8, abstol=1e-8)

println(sol.u[end][veh.id_r])
println(TVC.Utils.to_matrix(sol.u[end][veh.id_quat]))
println(sol.u[end][veh.id_ω])

plot(sol, vars=veh.id_v)
# plot(sol.t, veh.Acceleration.(sol.t))

#   Is optimiser initial constraint correct
#   ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡

traj=RocketTrajectoryParameters(r0=[22.366130464744867, -5.774597848558639, 28.23627501803703], v0=[3.8806161631580847, -2.910462122368562, -5.872678737716989], q0=[0.9969377999299179, -0.046919167783644675, -0.06255889037819291, -7.544599728408613e-20], ω0=[-0.3199687185611176, -0.4266249580814899, 2.3758547531861336e-17])

p = ODEParametersFeas(traj=traj, Control=(force=zeros(3), torque=zeros(3)), MotorIgnitionTime=solguidance.p[1])

x₀ = [traj.r0; traj.v0; traj.q0; traj.ω0][[p.veh.id_r; p.veh.id_v; p.veh.id_quat; p.veh.id_ω]];

#   Solve ODE
#   ≡≡≡≡≡≡≡≡≡≡≡

tspan = (traj.t0, p.MotorIgnitionTime)

prob = ODEProblem(f!, x₀, tspan, p)

soldiff = DifferentialEquations.solve(prob, reltol=1e-8, abstol=1e-8)

println(norm(soldiff(p.MotorIgnitionTime) - (solguidance.xd[1:13, 1])))

# traj = RocketTrajectoryParameters(q0=ones(4) / 2, ω0=zeros(3))
# t1 = ForwardDiff.derivative(t -> TVC.Utils.quatL(traj.q0) * TVC.Utils.wexp(traj.ω0 * t), 0.5)

# t2 = TVC.Utils.quatL(traj.q0) * [-sin(norm(traj.ω0 * t) / 2) * norm(traj.ω0) / 2; (traj.ω0 * t) / norm(traj.ω0 * t) * cos(norm(traj.ω0 * t) / 2) * norm(traj.ω0) / 2] # assume p[veh.id_tcoast] ≥ 0?

# norm(t1 - t2)




########## Test solving SCP Solution from previous solution
using TVC, SCPToolbox

veh = RocketParameters()
atmos = Atmosphere()
traj = RocketTrajectoryParameters();

mdl = RocketProblem(veh, atmos, traj)

# servoΔt = 0.02 # Servo step rate
x_0 = [traj.r0; traj.v0; traj.q0; traj.ω0][[veh.id_r; veh.id_v; veh.id_quat; veh.id_ω]];

#   Solve Optimal Trajectory
#   ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡

sol = solveProblem(mdl);

tₘ = 3.4239727124217527 
t0 = traj.t0 # 0
time = (tₘ - t0) / (veh.BurnTime - t0)

x = sample(solution.xc, time)

r = x[veh.id_r]
v = x[veh.id_v]
quat = x[veh.id_quat]
ω = x[veh.id_ω]
T = x[veh.id_T]
Ṫ = x[veh.id_Ṫ]

traj2 = RocketTrajectoryParameters(r0=r, v0=v, q0=quat, ω0=ω, T0=T, Ṫ0=Ṫ, t0=tₘ, MotorFired=true, PreviousTrajectoryState=sol.xc, PreviousTrajectoryInput=sol.uc, PreviousTrajectoryCurrentTime=time, UsePreviousTrajectory=true)

mdl2 = RocketProblem(veh, atmos, traj2)

sol2 = solveProblem(mdl2)

using Plots
plotly()

Plots.plot(traj.t0 .+ (mdl.veh.BurnTime - traj.t0) * sol.td, sol.xd[mdl.veh.id_r, :]')
Plots.plot!(traj2.t0 .+ (mdl2.veh.BurnTime - traj2.t0) * sol2.td, sol2.xd[mdl.veh.id_r, :]')

N = max(floor(Int, (mdl2.veh.BurnTime - mdl2.traj.t0) / 0.4), 0) + 2
SampleTimes = collect(range(traj2.PreviousTrajectoryCurrentTime, 1, N))
x = mapreduce(t -> sample(traj2.PreviousTrajectoryState, t), hcat, SampleTimes)
u = mapreduce(t -> sample(traj2.PreviousTrajectoryInput, t), hcat, SampleTimes)

Plots.plot!(traj2.t0 .+ (mdl2.veh.BurnTime - traj2.t0) * sol2.td, x[mdl.veh.id_r, :]')



tₘ = 3.4239727124217527 
t0 = traj.t0 # 0
t0 = tₘ - .25
time = (tₘ - t0) / (veh.BurnTime - t0)
time = 0.905707555921009

r = [32.500210374646585, -13.366532195806492, 0.02490321322506425]
v = [0.28295186872507366, -0.21110315458647022, -0.3071411144887891]
q = [0.9928946098146643, -9.906016038958462e-5, -0.00014485698746154275, -0.11896228659441592]
ω = [-0.0006710074693798252, -0.00045055087177835727, 0.14079149141098815]
T = [0.0031893090891592616, 0.0002654782354665, 0.9762048207819294]
Ṫ = [-0.024727024395630348, 0.02357795533709619, -0.03023221847387955]

traj3 = RocketTrajectoryParameters(r0=r, v0=v, q0=q, ω0=ω, T0=T, Ṫ0=Ṫ, t0=tₘ, MotorFired=true, PreviousTrajectoryState=sol.xc, PreviousTrajectoryInput=sol.uc, PreviousTrajectoryCurrentTime=time, UsePreviousTrajectory=true)


traj3 = RocketTrajectoryParameters(r0=r, v0=v, q0=q, ω0=ω, T0=T, Ṫ0=Ṫ, t0=tₘ, MotorFired=true)

mdl3 = RocketProblem(veh, atmos, traj3)

sol3, hist = solveProblem(mdl3);

N = max(floor(Int, (mdl3.veh.BurnTime - mdl3.traj.t0) / 0.4), 0) + 2
SampleTimes = collect(range(traj2.PreviousTrajectoryCurrentTime, 1, N))
x = mapreduce(t -> sample(traj2.PreviousTrajectoryState, t), hcat, SampleTimes)
u = mapreduce(t -> sample(traj2.PreviousTrajectoryInput, t), hcat, SampleTimes)

Plots.plot(traj3.t0 .+ (mdl3.veh.BurnTime - traj3.t0) * sol3.td, x[mdl.veh.id_r, :]')




##############

tₘ = 0.
t0 = traj.t0 # 0

r = [20.99451621999152, -4.745887164993641, 29.693564875235637]
v = [3.95530484742001, -2.9664786355650072, -2.4512205882479527]
q = [0.9999137738116913, -0.007879115279663097, -0.010505487039550798, 2.7982339831841635e-18]
ω = [-0.12744342006572731, -0.16992456008763648, 5.062291982155317e-16]
T = [0.0, 0.0, 1.0]
Ṫ = [0.0, 0.0, 0.0]

traj4 = RocketTrajectoryParameters(r0=r, v0=v, q0=q, ω0=ω, T0=T, Ṫ0=Ṫ, t0=tₘ, PreviousTrajectoryState=sol.xc, PreviousTrajectoryInput=sol.uc, PreviousTrajectoryCurrentTime=time, UsePreviousTrajectory=true)


traj4 = RocketTrajectoryParameters(r0=r, v0=v, q0=q, ω0=ω, T0=T, Ṫ0=Ṫ, t0=tₘ)

mdl4 = RocketProblem(veh, atmos, traj4)

sol4 = solveProblem(mdl4);

model = hist.subproblems[1].prg.mdl
grb = unsafe_backend(model)
# Gurobi.column(grb, x[1])
num_constraints = 708 #number of constraints? # Number of rows in constraint matrix
basis = Vector{Cint}(undef, num_constraints)
GRBgetBasisHead(grb, basis)
GRBgetconstrs(grb, basis)
GRBgetcoeff

GRBgetgenconstrMin(grb, )

numnzP = Ref{Cint}()
row = Cint(0) #Cint(Gurobi._info(model, c).row - 1)
ret = GRBgetconstrs(grb, numnzP, C_NULL, C_NULL, C_NULL, row, 1)
Gurobi._check_ret(model, ret)
cbeg = Array{Cint}(undef, 2)
cind = Array{Cint}(undef, numnzP[])
cval = Array{Cdouble}(undef, numnzP[])
ret = GRBgetconstrs(grb, numnzP, cbeg, cind, cval, row, 1)
Gurobi._check_ret(model, ret)

list_of_constraint_types(model)
all_constraints(model, list_of_constraint_types(model)[1]...)