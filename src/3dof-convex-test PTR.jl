using JuMP, ECOS, MathOptInterface, Utils, Plots, LinearAlgebra
using DifferentialEquations
include("./Rocket_Acceleration.jl")
include("./Quaternions.jl")

t_coast = 3.45
dt = 0.1 # see LanderSolid.m, dt needs to be low. We need many degrees of freedom on u

N = floor(Int, t_coast / dt) + 1 
dt = t_coast / (N - 1)

wvc = 1e3 # wtr needs to be small, otherwise we get suboptimal solution, probably a local optima?

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
    @objective(model, Min, t + wtr*(trapz(ηx, time)+trapz(ηu, time)) + wtr * ηp + wvc*(trapz(P, time))  )

    # @objective(model, Min, trapz(dot.(eachcol(r), eachcol(r)), collect(0:dt:t_coast)) ) # norm not supported in this fashion? cba to work out how

    @constraint(model, p >= 0)

    @constraint(model, r[:, 1] .== r_0 + p_prev * v_0 + p_prev^2 /2 * g + (v_0 + g * p_prev) * (p - p_prev))
    @constraint(model, v[:, 1] .== v_0 + p_prev * g + g * (p - p_prev))

    @constraint(model, r[3, N] == r_N[3])

    @constraint(model, r[3, :] .>= 0)

    # u1(k) = u[:, k] / norm(u[:, k])
    # u2(k) = u[:, k + 1] / norm(u[:, k + 1])

    # if slerpBool
    #     @constraint(model, [k in 1:N-1], v[:, k + 1] .== v[:, k] + g * dt + dt / θ^2 * ((T1 - T2) * (u1 - u2) - θ * cot(θ) * (T1 * u1 + T2 * u2) + θ * csc(θ) *  (T2 * u1 + T1 * u2)) + virtual[4:6, k])
    
    #     @constraint(model, [k in 1:N-1], r[:, k + 1] .== r[:, k] + v[:, k] * dt + g * dt^2 / 2 + dt^2 / θ^2 * (2 * T1 * u1 - T2 * (u1 + u2)) + dt^2 / θ^3 * cot(θ) * (-T1 * (-2 + θ^2) * u1 + 2 * T1 * u2 - 2 * T2 * (u1 + u2)) + dt^2 / θ^3 * csc(θ) * (-2 * T1 * u1 + T1 * (-2 + θ^2) * u2 + 2 * T2 * (u1 + u2)) + virtual[1:3, k])
    # else
        @constraint(model, [k in 1:N-1], v[:, k + 1] .== v[:, k] + (g + u[:, k] + g + u[:, k + 1]) / 2 * dt + virtual[4:6, k])
        @constraint(model, [k in 1:N-1], r[:, k + 1] .== r[:, k] + dt / 2 * (v[:, k] + v[:, k + 1]) + dt^2 / 12 * (u[:, k] - u[:, k + 1]) + virtual[1:3, k])
    # end
    
    @constraint(model, [k in 1:N], [Acceleration(dt * (k - 1)); u[:, k]] in SecondOrderCone()) # @constraint(model, [k in 1:N], norm(u[:, k]) <= 10)

    # Angular velocity constraint
    # Non convex constraint is u(k + 1) dot u(k) >= cos(w_max) * |u(k)| * |u(k + 1)|

    w_max = cos(deg2rad(45))

    @constraint(model, [[0; 0; 1]' * u[:, 2] / (w_max * norm([0; 0; 1])); u[:, 2]] in SecondOrderCone())

    @constraint(model, [[0; 0; 1]' * u[:, N - 1] / (w_max * norm([0; 0; 1])); u[:, N - 1]] in SecondOrderCone())

    f0(k) = u_prev[:, k + 1]' * u_prev[:, k] - w_max * (norm(u_prev[:, k + 1]) * norm(u_prev[:, k]))
    f_k(k) = u_prev[:, k + 1]' - w_max * norm(u_prev[:, k + 1]) * u_prev[:, k]' / norm(u_prev[:, k])
    f_k1(k) = u_prev[:, k]' - w_max * norm(u_prev[:, k]) * u_prev[:, k + 1]' / norm(u_prev[:, k + 1])

    @constraint(model, [k in 2:N-1], f0(k) + f_k(k) * (u[:, k] - u_prev[:, k]) + f_k1(k)  * (u[:, k + 1] - u_prev[:, k + 1]) >= 0)

    @constraint(model, delta_x .== [r; v] .- x_prev) 
    @constraint(model, delta_u .== u .- u_prev) 
    @constraint(model, delta_p .== p - p_prev) 

    return model
end

println("Augmented cost               Cost                  Δu                    Coast time (s)             Equality Slack")

for i = 1:10
    local wtr = i < 5 ? 1e-3 : 1e1 # if we don't do this u[:, 3] bounces around in order to break ω_max constraint, and if wtr is too high we don't get to a good solution.
    local model = build(x, u, dt, p, t_burn, wvc, wtr)
    optimize!(model)
    # @show objective_value(model)
    # @show norm(JuMP.value.(model.obj_dict[:P]), Inf)
    # @show norm(JuMP.value.(model.obj_dict[:ηx]), Inf)

    local r = JuMP.value.(model.obj_dict[:r])
    local v = JuMP.value.(model.obj_dict[:v])  
    global x = [r; v]

    global u_prev = u
    global u = JuMP.value.(model.obj_dict[:u])

    global p = JuMP.value.(model.obj_dict[:p])

    local virtual = JuMP.value.(model.obj_dict[:virtual])

    println("$(objective_value(model))    $(norm(v[:, N] .- v_N))    $(norm(JuMP.value.(model.obj_dict[:ηu]), Inf))             $p                  $(norm(virtual))")
end

# model = build(x, u, dt, p, t_burn, wvc, wtr)
# optimize!(model)
# @show objective_value(model)

# r = JuMP.value.(model.obj_dict[:r])
# v = JuMP.value.(model.obj_dict[:v])  
# u = JuMP.value.(model.obj_dict[:u])
# p = JuMP.value.(model.obj_dict[:p])

# x = [r; v]

t = 0:dt:t_burn;

plotlyjs()

plot(t, norm.(eachcol(u)))
plot(t, ThrottleLevel(norm.(eachcol(u)), Acceleration(t) ) )

u_normalised = deepcopy(u)

for k = 1:N
    if norm(u_normalised[:, k]) <= 2
        u_normalised[:, k] = [0; 0; 1]
    end
end
ω_max = zeros(N - 1)
u_normalised =  u_normalised ./ norm.(eachcol(u_normalised))'
# test = u ./ norm.(eachcol(u))';
# acosd(test[:, 2]' * test[:, 3])
for k = 1:N-1
    ω_max[k] = acosd(u_normalised[:, k]' * u_normalised[:, k + 1])
end
plot(t[2:N], ω_max)

# Continuous version of u
function uc(u, u_normalised, dt, t, slerpBool, realAcceleration) 
    if t < 0 || t >= t_burn 
        return [0; 0; 0]
    else 
        tmp = Int(t ÷ dt) #floor(Int, t / dt)
        frac = (t % dt) / dt

        if slerpBool
            u_norm = slerp(u_normalised[:, tmp + 1], u_normalised[:, tmp + 2], frac)
            normu = norm(u[:, tmp + 1]) + frac * (norm(u[:, tmp + 2]) - norm(u[:, tmp + 1]))
            if realAcceleration
                return u_norm * Acceleration(t)
            else
                return u_norm * normu
            end
        else
            return u[:, tmp + 1] + frac * (u[:, tmp + 2] - u[:, tmp + 1])
        end
    end
end

function f!(dx, x, p, t)
    dx[1:3] = x[4:6]
    Thrust = uc(u, u_normalised, dt, t, true, true) 
    dx[4:6] = g + Thrust
end

prob = ODEProblem(f!, x[:, 1], (0.0, t_burn))
sol = DifferentialEquations.solve(prob, reltol=1e-12, abstol=1e-12)

xc = sol(0:dt:t_burn).u
xmatrix = transpose(reduce(hcat, xc))
# plot(xmatrix[:, 1], xmatrix[:, 2], xmatrix[:, 3])
# # plot(sol)

println(norm(xmatrix - x', Inf))
println(norm(xmatrix - x'))

println("Cost: $(norm(xmatrix[end, 4:6] - v_N))")

# Same as when slerpBool = false.

# r = zeros(3, N)
# v = zeros(3, N)
# r[:, 1] = x[1:3, 1]
# v[:, 1] = x[4:6, 1]
# for k = 1:N-1
#     v[:, k + 1] = v[:, k] + (g + u[:, k] + g + u[:, k + 1]) / 2 * dt
#     r[:, k + 1] = r[:, k] + dt / 2 * (v[:, k] + v[:, k + 1]) + dt^2 / 12 * (u[:, k] - u[:, k + 1])
#     # r[: , 2] = x[1:3, 1] + x[4:6] * t + t^2/2 * (g + u[:, 1]) + t^3 / (6 * dt) * (u[:, 2] - u[:, 1]) # dt is time between u samples.
# end
# x_exact = [r; v]

# Same as when slerpBool = true.

r = zeros(3, N)
v = zeros(3, N)

r[:, 1] = x[1:3, 1]
v[:, 1] = x[4:6, 1]

Δt = t_coast / (N - 1)
for k = 1:N-1
    θ = acos(u_normalised[:, k]' * u_normalised[:, k + 1]) 
    # local t = 0 + (k - 1) * Δt
    # frac = (t % dt) / dt
    u1 = u_normalised[:, k]
    u2 = u_normalised[:, k + 1]
    T1 = norm(u[:, k])
    T2 = norm(u[:, k + 1])
    
    # https://www.wolframalpha.com/input/?i=integrate+%28sin%28q+*+%281+-+t%2Fh%29%29+*+u_1+%2B+sin%28q+*+t+%2F+h%29+*+u_2%29+*+%28A+%2B+t+*+B%29, the  / sin(theta), g * t terms are missing h = dt, B = (T_2 - T_1) / dt, q = theta
    # v[:, k + 1] = v[:, k] + g * Δt + 1 / (θ * sin(θ)) * (T1 * dt * (cos(θ * (1 - frac)) * u1 - cos(θ * frac) * u2) + (T2 - T1) * dt / θ * (sin(θ * (1 - frac)) * u1 + sin(θ * frac) * u2) + (T2 - T1) * Δt * (cos(θ * (1 - frac)) * u1 - cos(θ * frac) * u2)) # Δt is timestep for sim, dt is for control interpolation, no need for frac here.

    # v[:, k + 1] = v[:, k] + g * dt + 1 / (θ * sin(θ)) * (T1 * dt * (u1 - cos(θ) * u2) + (T2 - T1) * dt / θ * sin(θ) * u2 + (T2 - T1) * dt * (u1 - cos(θ) * u2)) # assume frac = 1

    # dt / (θ * sin(θ)) * (T1 * (u1 - cos(θ) * u2) + (T2 - T1) * (sin(θ) / θ * u2 + u1 - cos(θ) * u2)) # further simplified


    # https://www.wolframalpha.com/input/?i=integrate+%28sin%28q+*+%281+-+t%2Fh%29%29+*+u_1+%2B+sin%28q+*+t+%2F+h%29+*+u_2%29+*+%28A+%2B+t+*+B%29+from+t+%3D0+to+t+%3D+h
    # v[:, k + 1] = v[:, k] + g * dt + dt / (θ * sin(θ)) * ( u1 * (T1 * θ * (1 - cos(θ)) + (T2 - T1) * (θ - sin(θ))) + u2 * (-θ * cos(θ) * T2 + T1 * θ + (T2 - T1) * sin(θ)) )


    # https://www.wolframalpha.com/input/?i=integrate+%28sin%28q+*+%281+-+t%2Fh%29%29+*+u_1+%2B+sin%28q+*+%28t%2Fh%29%29+*+u_2%29+*+%281+-+t%2Fh%29+*+T_1+%2F+sin%28q%29+from+t+%3D0+to+t+%3D+h
    # https://www.wolframalpha.com/input/?i=integrate+%28sin%28q+*+%281+-+t%2Fh%29%29+*+u_1+%2B+sin%28q+*+%28t%2Fh%29%29+*+u_2%29+*+%28t%2Fh%29+*+T_2+%2F+sin%28q%29+from+t+%3D0+to+t+%3D+h
    # v[:, k + 1] = v[:, k] + g * dt - dt / θ^2 * T1 * (u1 * (θ * cot(θ) - 1) + u2 * (1 - θ * csc(θ))) - dt / θ^2 * T2 * u2 * ( θ * cot(θ) - 1)
    # v[:, k + 1] = v[:, k] + g * dt - dt / θ^2 * T1 * u1 * (- θ * tan(θ/2) + θ * csc(θ) - 1) - dt / θ^2 * T2 * u2 * ( θ * cot(θ) - 1)
    

    # Check mathematica notebook
    v[:, k + 1] = v[:, k] + g * dt + dt / θ^2 * ((T1 - T2) * (u1 - u2) - θ * cot(θ) * (T1 * u1 + T2 * u2) + θ * csc(θ) *  (T2 * u1 + T1 * u2))
    
    r[:, k + 1] = r[:, k] + v[:, k] * dt + g * dt^2 / 2 + dt^2 / θ^2 * (2 * T1 * u1 - T2 * (u1 + u2)) + dt^2 / θ^3 * cot(θ) * (-T1 * (-2 + θ^2) * u1 + 2 * T1 * u2 - 2 * T2 * (u1 + u2)) + dt^2 / θ^3 * csc(θ) * (-2 * T1 * u1 + T1 * (-2 + θ^2) * u2 + 2 * T2 * (u1 + u2))
end
x_slerp_exact = [r; v];

println(norm(x_slerp_exact - x, Inf))
println(norm(x_slerp_exact - x))

println("Cost with slerp: $(norm(x_slerp_exact[4:6, end] - v_N))")

# println(norm(x_exact - x, Inf))
# println(norm(x_exact - x))

# r = zeros(3, N)
# v = zeros(3, N) 

# u = -g + tmp

# r[:, 1] = r_0
# v[:, 1] = v_0

# for k in 1:N-1
#     v[:, k + 1] = v[:, k] + (g + u) * dt
#     r[:, k + 1] = r[:, k] + dt / 2 * (v[:, k] + v[:, k + 1])
# end