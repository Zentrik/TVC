using JuMP, ECOS

N = 11
t_coast = 3.45
dt = t_coast / (N - 1)

g = [0; 0; -9.80655]

r_0 = [20; -4; 30];	# position vector, m
v_0 = [4; -3; 0];		# velocity vector, m/s
r_N =[0; 0; 0];				# terminal position, m
v_N =[0; 0; 0];	# terminal velocity, m

model = Model(ECOS.Optimizer) 

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
@constraint(model, [k in 1:N-1], r[:, k + 1] .== r[:, k] + dt / 2 * (v[:, k] + v[:, k + 1]) + dt^2 / 12 * (u[:, k + 1] - u[:, k]) )

@constraint(model, [k in 1:N], [10; u[:, k]] in SecondOrderCone()) # @constraint(model, [k in 1:N], norm(u[:, k]) <= 10)

#print(model)
optimize!(model)
@show objective_value(model)

u = JuMP.value.(u)
v  = JuMP.value.(v)

t = 0:dt:t_coast;
plot(t, norm.(eachcol(v)))
plot(t, norm.(eachcol(u)))

# r = zeros(3, N)
# v = zeros(3, N) 

# u = -g + tmp

# r[:, 1] = r_0
# v[:, 1] = v_0

# for k in 1:N-1
#     v[:, k + 1] = v[:, k] + (g + u) * dt
#     r[:, k + 1] = r[:, k] + dt / 2 * (v[:, k] + v[:, k + 1])
# end