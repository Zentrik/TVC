using JuMP, ECOS, MathOptInterface, Utils, Plots, LinearAlgebra
include("./Rocket_Acceleration.jl")

t_coast = 3.45
dt = 0.1 # see LanderSolid.m, dt needs to be low. We need many degrees of freedom on u

N = floor(Int, t_coast / dt) + 1 
dt = t_coast / (N - 1)

wvc, wtr = 1e3, 1e-4 # wtr needs to be small, otherwise we get suboptimal solution, probably a local optima?

g = [0; 0; -9.80655]

r_0 = [20; -4; 30];	# position vector, m
v_0 = [4; -3; 0];		# velocity vector, m/s
r_N =[0; 0; 0];				# terminal position, m
v_N =[0; 0; 0];	# terminal velocity, m

t_coast = 0.0
t_burn = 3.45
tmp = [0; 0; - 2 * (r_0[3] + v_0[3] * t_burn) / t_burn^2] # [0.0, 0.0, -5.040957781978575]

x = zeros(6, N)
u = zeros(3, N)
p = t_coast

for k = 1:N
    local t = (k - 1) / (N - 1) * t_burn
    x[:, k] = [r_0 + t * v_0 + t^2 * tmp / 2; v_0 + tmp * t]
    u[:, k] = -g + tmp # [0.0, 0.0, 4.765592218021425]
end   

#print(model)
function build(x_prev, u_prev, dt, p_prev, t_burn, wvc, wtr)

    model = Model(optimizer_with_attributes(ECOS.Optimizer, "verbose" => 0)); 

    @variable(model, r[1:3, 1:N])
    @variable(model, v[1:3, 1:N])

    @variable(model, u[1:3, 1:N])

    @variable(model, p)

    @variable(model, virtual[1:6, 1:N])

    @variable(model, delta_x[1:6, 1:N])
    @variable(model, delta_u[1:3, 1:N])
    @variable(model, delta_p)

    t = @variable(model)
    @constraint(model, [t; v[:, N] .- v_N] in SecondOrderCone()) # norm objective

    # qtr = Inf
    @variable(model, ηx[1:N])
    @constraint(model, [k in 1:N], [ηx[k]; delta_x[:, k]] in MathOptInterface.NormInfinityCone(7))

    @variable(model, ηu[1:N])
    @constraint(model, [k in 1:N], [ηu[k]; delta_u[:, k]] in MathOptInterface.NormInfinityCone(4))

    @variable(model, ηp)
    @constraint(model, [ηp; delta_p] in MathOptInterface.NormInfinityCone(2))

    @variable(model, P[1:N])
    @constraint(model, [k in 1:N], [P[k]; virtual[:, k]] in MathOptInterface.NormOneCone(7))

    time = collect(0:dt:t_burn)
    @objective(model, Min, t + wtr*(trapz(ηx, time)+trapz(ηu, time)) + 0.1 * ηp + wvc*(trapz(P, time))  )

    # @objective(model, Min, trapz(dot.(eachcol(r), eachcol(r)), collect(0:dt:t_coast)) ) # norm not supported in this fashion? cba to work out how

    @constraint(model, p >= 0)

    @constraint(model, r[:, 1] .== r_0 + p_prev * v_0 + p_prev^2 /2 * g + (v_0 + g * p_prev) * (p - p_prev))
    @constraint(model, v[:, 1] .== v_0 + p_prev * g + g * (p - p_prev))

    @constraint(model, r[3, N] == r_N[3])

    @constraint(model, r[3, :] .>= 0)

    @constraint(model, [k in 1:N-1], v[:, k + 1] .== v[:, k] + (g + u[:, k] + g + u[:, k + 1]) / 2 * dt + virtual[4:6, k])
    @constraint(model, [k in 1:N-1], r[:, k + 1] .== r[:, k] + dt / 2 * (v[:, k] + v[:, k + 1]) + dt^2 / 12 * (u[:, k + 1] - u[:, k]) + virtual[1:3, k])

    @constraint(model, [k in 1:N], [Acceleration(dt * (k - 1)); u[:, k]] in SecondOrderCone()) # @constraint(model, [k in 1:N], norm(u[:, k]) <= 10)

    @constraint(model, delta_x .== [r; v] .- x_prev) 
    @constraint(model, delta_u .== u .- u_prev) 
    @constraint(model, delta_p .== p - p_prev) 

    return model
end

for i = 1:0
    local model = build(x, u, dt, p, t_burn, wvc, wtr)
    optimize!(model)
    @show objective_value(model)
    @show norm(JuMP.value.(model.obj_dict[:P]), Inf)
    @show norm(JuMP.value.(model.obj_dict[:ηx]), Inf)

    local r = JuMP.value.(model.obj_dict[:r])
    local v = JuMP.value.(model.obj_dict[:v])  
    global x = [r; v]

    global u = JuMP.value.(model.obj_dict[:u])

    global p = JuMP.value.(model.obj_dict[:p])
end

model = build(x, u, dt, p, t_burn, wvc, wtr)
optimize!(model)
@show objective_value(model)

r = JuMP.value.(model.obj_dict[:r])
v = JuMP.value.(model.obj_dict[:v])  
u = JuMP.value.(model.obj_dict[:u])
p = JuMP.value.(model.obj_dict[:p])

x = [r; v]

t = 0:dt:t_burn;

plotlyjs()

plot(t, norm.(eachcol(u)))
plot(t, ThrottleLevel(norm.(eachcol(u)), Acceleration(t) ) )

print(p)

# r = zeros(3, N)
# v = zeros(3, N) 

# u = -g + tmp

# r[:, 1] = r_0
# v[:, 1] = v_0

# for k in 1:N-1
#     v[:, k + 1] = v[:, k] + (g + u) * dt
#     r[:, k + 1] = r[:, k] + dt / 2 * (v[:, k] + v[:, k + 1])
# end