using DifferentialEquations
include("./Rocket_Acceleration.jl")

g = [0; 0; -9.80655];

function rocket!(x_dot, x, p ,t)
    x_dot[1:3] = x[4:6]
    x_dot[4:6] = g + [0; 0; Acceleration(t)]
end

tspan = (0.0, 100.0) # should be long enough for rocket to hit ground

x0 = [20.0; -4.0; 30.0; 4.0; -3.0; 0.0]
prob = ODEProblem(rocket!,x0,tspan) # prob = ODEProblem(rocket!,[20; -4; 30;4; -3; 0], tspan, u0)

condition(x,t,integrator) = x[3] # when zero halt integration
affect!(integrator) = terminate!(integrator)
cb = ContinuousCallback(condition,affect!)
sol = DifferentialEquations.solve(prob,Tsit5(),callback=cb)

x = zeros(pbm.nx, N)
u = zeros(pbm.nu, N)

for k = 1:N
    t = (k - 1) / (N - 1) * tf
    x[:, k] = [r_0 + t * v_0 + t^2 * tmp / 2; v_0 + tmp * t]
    u[:, k] = u
end