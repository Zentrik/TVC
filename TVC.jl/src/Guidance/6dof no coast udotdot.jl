# # TODO: Make quaternion discretisation more accurate. Without constraint on quaternion norm, the norm of the quaternion strays from 1, I think this is due to the Jacobians of dynamics being evaluated about a reference trajectory that has non unit quaternions.

# using SCPToolbox

# using LinearAlgebra
# # using ECOS

# using ..Utils
# import ..Utils: rotate, skew
# using ..Guidance

# export define_problem_no_coast!

# # Per Successive Convexification for Mars 6-DoF Powered Descent Landing Guidance, 2017. Set control to second derivative of thrust vector
# # This allows me to add a constraint on the tvc gimbal rate and it adds more degrees of freedom on the control/ allows for more complex control for a given N (no. of discretisation steps) improving cost.

# function define_problem_no_coast!(pbm::TrajectoryProblem)::Nothing
#     set_dims!(pbm)
#     set_scale!(pbm)
#     set_cost!(pbm)
#     set_dynamics!(pbm)
#     set_integration_action(pbm)
#     set_convex_constraints!(pbm)
#     set_bcs!(pbm)

#     set_guess!(pbm)

#     return nothing
# end

# function set_dims!(pbm::TrajectoryProblem)::Nothing

#     problem_set_dims!(pbm, 19, 4, 1)

#     return nothing
# end

# function set_scale!(pbm::TrajectoryProblem)::Nothing #VERY IMPORTANT
#     advise! = problem_advise_scale!

#     # States
#     advise!(pbm, :state, 1, (-1000.0, 1000.0))
#     advise!(pbm, :state, 2, (-1000.0, 1000.0))
#     advise!(pbm, :state, 3, (0.0, 1000.0))
#     advise!(pbm, :state, 4, (-100.0, 100.0))
#     advise!(pbm, :state, 5, (-100.0, 100.0))
#     advise!(pbm, :state, 6, (-100.0, 100.0))

#     advise!(pbm, :state, 7, (-1.0, 1.0))
#     advise!(pbm, :state, 8, (-1.0, 1.0))
#     advise!(pbm, :state, 9, (-1.0, 1.0))
#     advise!(pbm, :state, 10, (-1.0, 1.0))
#     advise!(pbm, :state, 11, (-10.0, 10.0))
#     advise!(pbm, :state, 12, (-10.0, 10.0))
#     advise!(pbm, :state, 13, (-10.0, 00.0))

#     advise!(pbm, :state, 14, (-1.0, 1.0))
#     advise!(pbm, :state, 15, (-1.0, 1.0))
#     advise!(pbm, :state, 16, (-1.0, 1.0))
#     advise!(pbm, :state, 17, (-float(π), float(π)))
#     advise!(pbm, :state, 18, (-float(π), float(π)))
#     advise!(pbm, :state, 19, (-float(π), float(π)))

#     # Inputs
#     advise!(pbm, :input, 1, (-deg2rad(10), deg2rad(10)))
#     advise!(pbm, :input, 2, (-deg2rad(10), deg2rad(10)))
#     advise!(pbm, :input, 3, (-deg2rad(10), deg2rad(10)))
#     advise!(pbm, :input, 4, (-1.0, 1.0))

#     # Parameters
#     advise!(pbm, :parameter, 1, (0.0, 1.0))

#     return nothing
# end

# function straightline_interpolate(v0, vf, N::Int)

#     # Initialize
#     nv = length(v0)
#     v = zeros(nv, N)

#     for k = 1:N
#         mix = (k - 1) / (N - 1)

#         v[:, k] = v0 + mix * (vf - v0)
#     end

#     return v
# end

# function set_guess!(pbm::TrajectoryProblem)::Nothing
#     problem_set_guess!(pbm, (N, pbm) -> begin
#         veh = pbm.mdl.veh
#         traj = pbm.mdl.traj
#         atmos = pbm.mdl.atmos

#         p = [0.0] # ignite immediately.

#         # motorTimeRemaining = veh.BurnTime - traj.t0 # how much motor time remaining

#         # dt = motorTimeRemaining / (N - 1) # for convex constraints

#         x = zeros(pbm.nx, N)
#         u = zeros(pbm.nu, N)

#         x[veh.id_r, :] = straightline_interpolate(traj.r0, traj.rN, N)
#         x[veh.id_v, :] = straightline_interpolate(traj.v0, traj.vN, N)
#         x[veh.id_ω, :] = straightline_interpolate(traj.ω0, traj.ωN, N)

#         x[veh.id_T, :] = straightline_interpolate(traj.T0, traj.TN, N)  
#         x[veh.id_Ṫ, :] = straightline_interpolate(traj.Ṫ0, traj.ṪN, N) 
        
#         u = straightline_interpolate([0; 0; 0; 0], [0; 0; 0; 0], N)
        
#         for k = 1:N
#             mix = (k - 1) / (N - 1)
            
#             x[veh.id_quat, k] = slerp_quat(traj.q0, traj.qN, mix)
#         end

#         return x, u, p
#     end)

#     return nothing
# end
    
# function set_cost!(pbm::TrajectoryProblem)::Nothing
#     problem_set_terminal_cost!(
#         pbm, (x, p, pbm) -> dot(x[pbm.mdl.veh.id_v] - pbm.mdl.traj.vN, x[pbm.mdl.veh.id_v] - pbm.mdl.traj.vN)
#     )

#     return nothing
# end

# function set_dynamics!(pbm::TrajectoryProblem)::Nothing
#     function f_t(x, u, t, pbm)
#         veh = pbm.mdl.veh
#         atmos = pbm.mdl.atmos

#         r = x[veh.id_r]
#         v = x[veh.id_v]
#         quat = x[veh.id_quat]
#         ω = x[veh.id_ω]
#         T = x[veh.id_T]
#         Ṫ = x[veh.id_Ṫ]

#         return [
#         v; 
#         atmos.g(r[3]) + rotate(quat, T) * veh.Acceleration(t); 
#         1/2 * quatL(quat) * [0; ω]; 
#         veh.InertiaTensor \ (cross(veh.MomentArm(t), T * veh.Thrust(t)) + [0; 0; u[veh.id_roll]] - cross(ω, veh.InertiaTensor * ω)); 
#         Ṫ; 
#         u[veh.id_T̈]]
#     end

#     function A_t(x, u, t, veh)
#         r = x[veh.id_r]
#         v = x[veh.id_v]
#         quat = x[veh.id_quat]
#         ω = x[veh.id_ω]
#         T = x[veh.id_T]
#         Ṫ = x[veh.id_Ṫ]

#         A = zeros(pbm.nx, pbm.nx)
#         A[veh.id_r, veh.id_v] = I(3)
#         A[veh.id_v, veh.id_quat[1]] = veh.Acceleration(t) * 2 * (quat[1] * T + quat[2:4] × T)
#         A[veh.id_v, veh.id_quat[2:4]] = veh.Acceleration(t) * 2 * (quat[2:4]' * T * I(3) + quat[2:4] * T' - T * quat[2:4]' - quat[1] * skew(T))
#         A[veh.id_v, veh.id_T] = veh.Acceleration(t) * to_matrix(quat)
#         A[veh.id_quat, veh.id_quat] = 1/2 * quatR([0; ω])
#         A[veh.id_quat, veh.id_ω] = 1/2 * quatL(quat)[:, 2:4]
#         A[veh.id_ω, veh.id_ω] = - veh.InertiaTensor \ (skew(ω) * veh.InertiaTensor - skew(veh.InertiaTensor * ω))
#         A[veh.id_ω, veh.id_T] = veh.InertiaTensor \ skew(veh.MomentArm(t)) * veh.Thrust(t)
#         A[veh.id_T, veh.id_Ṫ] = I(3)

#         return A
#     end

#     function B_t(x, u, t, veh)
#         B = zeros(pbm.nx, pbm.nu)
#         B[veh.id_ω, veh.id_roll] = inv(veh.InertiaTensor)[:, 3]
#         B[veh.id_Ṫ, veh.id_T̈] = I(3)

#         return B
#     end

#     # x = rand(19, 1);
#     # u = rand(4, 1);
#     # Thrust = rand()
#     # count(f(x, u, Thrust) - f_t(x, u, Thrust) .> 1e-10)
#     # count(A(x, u, Thrust) - A_t(x, u, Thrust) .> 1e-10)
#     # count(B(x, u, Thrust) - B_t(x, u, Thrust) .> 1e-10)

#     # Dynamics
#     problem_set_dynamics!(
#         pbm,
#         # f
#         (t, k, x, u, p, pbm) -> f_t(x, u, motorTime(t, pbm.mdl), pbm) * (pbm.mdl.veh.BurnTime - pbm.mdl.traj.t0),
#         # df/dx
#         (t, k, x, u, p, pbm) -> A_t(x, u, motorTime(t, pbm.mdl), pbm.mdl.veh) * (pbm.mdl.veh.BurnTime - pbm.mdl.traj.t0),
#         # df/du
#         (t, k, x, u, p, pbm) -> B_t(x, u, motorTime(t, pbm.mdl), pbm.mdl.veh) * (pbm.mdl.veh.BurnTime  - pbm.mdl.traj.t0),
#         # df/dp
#         (t, k, x, u, p, pbm) -> zeros(pbm.nx, pbm.np) * (pbm.mdl.veh.BurnTime  - pbm.mdl.traj.t0)
#     )

#     return nothing
# end

# function set_integration_action(pbm::TrajectoryProblem)::Nothing

#     # Quaternion re-normalization on numerical integration step
#     problem_set_integration_action!(
#         pbm, pbm.mdl.veh.id_quat,
#         (q, pbm) -> begin
#             qn = q / norm(q)
#             return qn
#     end)

#     return nothing
# end

# function set_bcs!(pbm::TrajectoryProblem)::Nothing
#     # Boundary conditions
    
#     problem_set_bc!(
#         pbm, :ic, # Initial condition
#         (x, p, pbm) -> begin
#             veh = pbm.mdl.veh
#             traj = pbm.mdl.traj

#             x0 = zeros(pbm.nx)
#             x0[veh.id_r] = traj.r0
#             x0[veh.id_v] = traj.v0
#             x0[veh.id_quat] = traj.q0
#             x0[veh.id_ω] = traj.ω0
#             x0[veh.id_T] = traj.T0
#             x0[veh.id_Ṫ] = traj.Ṫ0

#             return x - x0
#         end,
#         (x, p, pbm) -> I(pbm.nx), # Jacobian wrt x
#         (x, p, pbm) -> zeros(pbm.nx, pbm.np) # Jacobian wrt p
#     )

#     problem_set_bc!(
#             pbm, :tc, # Terminal condition
#             (x, p, pbm) -> begin
#                 veh = pbm.mdl.veh
#                 traj = pbm.mdl.traj

#                 xf = zeros(9)
#                 xf[1] = traj.rN[3]
#                 xf[2:3] = traj.qN[2:3]
#                 xf[4:6] = traj.ωN
#                 xf[7:9] = traj.TN

#                 return x[vcat(veh.id_r[3], veh.id_quat[2:3], veh.id_ω, veh.id_T)] - xf
#             end,
#             (x, p, pbm) -> begin # Jacobian wrt x 
#                 veh = pbm.mdl.veh
#                 traj = pbm.mdl.traj

#                 J = zeros(9, pbm.nx)
#                 J[1, veh.id_r[3]] = 1
#                 J[2:3, veh.id_quat[2:3]] = I(2)
#                 J[4:6, veh.id_ω] = I(3)
#                 J[7:9, veh.id_T] = I(3)

#                 return J
#             end,
#             (x, p, pbm) -> zeros(9, pbm.np), # Jacobian wrt p
#         )
# end

# function set_convex_constraints!(pbm::TrajectoryProblem)::Nothing
#     # Convex State Constraints
#     problem_set_X!(
#         pbm, (t, k, x, p, pbm, ocp) -> begin
#             veh = pbm.mdl.veh
#             traj = pbm.mdl.traj

#             @add_constraint(
#                 ocp, NONPOS, "height >= 0", (x[veh.id_r[3]],), begin
#                 local height = arg[1]
#                 - height
#                 end)

#             @add_constraint(
#                 ocp, SOC, "Thrust Magnitude <= Max", (x[veh.id_T],), begin # we have say thrust <= 1, as we want it normalised
#                     local Thrust = arg[1]
#                     [1; Thrust]
#                 end)

#             @add_constraint(
#                 ocp, SOC, "Thrust Gimal angle < delta_max", (x[veh.id_T],), begin
#                     local Thrust = arg[1]
#                     [Thrust[3] / cos(pi * 5/180); Thrust]
#                 end)
            
#             @add_constraint( # shouldn't be necessary, but seems very very useful.
#                 ocp, SOC, "|quat| <= 1", (x[veh.id_quat],), begin
#                     local quat = arg[1]
#                     [1; quat]
#                 end)

#             # @add_constraint(
#             # ocp, SOC, "TVC angular velocity <= Max", (x[veh.id_Ṫ],), begin
#             #     local u = arg[1]
#             #     [deg2rad(5); u] # Angular velocity is r × v / ||r||², v = u, assume u is ⊥ r and ||r||² = 1, so angular velocity is v = u. 
#             # end)

#             @add_constraint(
#             ocp, SOC, "TVC angular velocity <= Max", (x[veh.id_Ṫ],), (u) -> [deg2rad(5); u]) # Angular velocity is r × v / ||r||², v = u, assume u is ⊥ r and ||r||² = 1, so angular velocity is v = u.

#         # @add_constraint(
#         #     ocp, SOC, "|w| < w_max", (x[11:13],), begin
#         #         local w = arg[1]
#         #         [pi / 2; w]
#         #     end)
#     end)

#     # Convex Input Constraints
#     problem_set_U!(
#         pbm, (t, k, u, p, pbm, ocp) -> begin
#             veh = pbm.mdl.veh
#             traj = pbm.mdl.traj

#             @add_constraint(
#                 ocp, L1, "|roll torque| < max", (u[4],), begin
#                     local u = arg[1]
#                     [0.1; u]
#                 end)

#             @add_constraint(
#                 ocp, SOC, "TVC Acceleration <= Max", (u[veh.id_T̈],), begin
#                     local u = arg[1]
#                     [deg2rad(10); u] # Angular Acceleration is r × a / ||r||², a = u, assume u is ⊥ r and ||r||² = 1, so angular acceleration is a = u. 
#                 end)
#     end) 
# end