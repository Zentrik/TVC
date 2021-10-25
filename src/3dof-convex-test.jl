using JuMP, ECOS, LinearAlgebra, Plots
include("./Rocket_Acceleration.jl")

t_burn = 3.45
dt = 0.1 # see LanderSolid.m, dt needs to be low. We need many degrees of freedom on u

N = floor(Int, t_burn / dt) + 1 
dt = t_burn / (N - 1)

g = [0; 0; -9.80655]

r0 = [20; -4; 30]	# position vector, m
v0 = [4; -3; 0] # velocity vector, m/s
r_N =[0; 0; 0];				# terminal position, m
v_N =[0; 0; 0];	# terminal velocity, m

list = 0:0.05:1.33
values = zeros(length(list))
for (i, t_coast) = enumerate(list)

    r_0 = [20; -4; 30] + v0 * t_coast + g * t_coast^2/2;	# position vector, m
    v_0 = [4; -3; 0] + g * t_coast;		# velocity vector, m/s

    model = Model(ECOS.Optimizer) 
    MOI.set(model, MOI.Silent(), true)

    @variable(model, r[1:3, 1:N])
    @variable(model, v[1:3, 1:N])
    @variable(model, u[1:3, 1:N])

    t = @variable(model)
    @constraint(model, [t; v[:, N] .- v_N] in SecondOrderCone()) # norm objective
    @objective(model, Min, t)

    # @objective(model, Min, trapz(dot.(eachcol(r), eachcol(r)), collect(0:dt:t_coast)) ) # norm not supported in this fashion? cba to work out how

    @constraint(model, r[:, 1] .== r_0)
    @constraint(model, v[:, 1] .== v_0)

    @constraint(model, r[3, N] == r_N[3])

    @constraint(model, r[3, :] .>= 0)

    @constraint(model, [k in 1:N-1], v[:, k + 1] .== v[:, k] + (g + u[:, k] + g + u[:, k + 1]) / 2 * dt )
    # v[:, k] + dt * g + dt / 2 * (u[:, k]
    @constraint(model, [k in 1:N-1], r[:, k + 1] .== r[:, k] + dt / 2 * (v[:, k] + v[:, k + 1]) + dt^2 / 12 * (u[:, k + 1] - u[:, k]) )
    # r[:, k] + dt  * v[:, k] + dt^2 / 2 * g + dt^2 / 3 * u[:, k] + dt^2 / 6 * u[:, k+1] 
    @constraint(model, [k in 1:N], [Acceleration(dt * (k - 1)); u[:, k]] in SecondOrderCone()) # @constraint(model, [k in 1:N], norm(u[:, k]) <= 10)

    #print(model)
    optimize!(model)

    values[i] = objective_value(model)
end

print(transpose(values))
plot(list, values)

# r = zeros(3, N)
# v = zeros(3, N) 

# u = -g + tmp

# r[:, 1] = r_0
# v[:, 1] = v_0

# for k in 1:N-1
#     v[:, k + 1] = v[:, k] + (g + u) * dt
#     r[:, k + 1] = r[:, k] + dt / 2 * (v[:, k] + v[:, k + 1])
# end