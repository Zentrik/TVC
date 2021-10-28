module scp_new_problem
include("./3dof t_burn more than 3.45.jl")
end # module

using .scp_new_problem

solution = scp_new_problem.solve(:ptr); # Remember J is augmented cost function
# print(sample(solution.xc, 1)[4:6]' * sample(solution.xc, 1)[4:6])

include("./Rocket_Acceleration.jl")
using Utils, LinearAlgebra, Plots


ctres = 1000
N = size(solution.xd, 2)
xct = hcat([sample(solution.xc, t) for t in LinRange(0, 1, ctres)]...)
vct = hcat([sample(solution.uc, t) for t in LinRange(0, 1, ctres)]...)

plotlyjs()

t_burn = length(solution.p) >= 2 ? solution.p[2] : 3.45

# Plots.plot(LinRange(0, 1, ctres) .* t_burn, xct[1:3, :]', title = "Position")

# Plots.plot(LinRange(0, 1, ctres) .* t_burn, xct[4:6, :]', title = "Velocity")

# Plots.plot(LinRange(0, 1, ctres) .* t_burn, vct')
plot(LinRange(0, 1, ctres) .* t_burn, Acceleration(LinRange(0, 1, ctres) .* t_burn))
plot!(LinRange(0, 1, ctres) .* t_burn, norm.(eachcol(vct)), title = "Thrust Magnitude") # norm constraint is broken intersample ...
# plot(solution.td .* t_coast, norm.(eachcol(solution.ud)) ./ 10.0, title = "Throttle Level")

# T = LinRange(0, 1, ctres) .* t_burn
# Magnitude = norm.(eachcol(vct));
# Max = Acceleration(T)

# # Continouous 
# plot(T, ThrottleLevel(Magnitude, Max), title = "Throttle Level")

# # Discrete
# plot(solution.td .* t_burn, ThrottleLevel(norm.(eachcol(solution.ud)), Acceleration(solution.td .* t_burn)), title = "Throttle Level")


println("Coast time (s): ", solution.p[1])
println("Burn time (s): ", t_burn)

v_N =[0; 0; 0];
println("Impact Velocity Magnitude (m/s): ", solution.cost^0.5)


# using Solvers
# fieldnames(Solvers.SCPSolution)

# outfile = "X.txt"
# open(outfile, "w") do f
#   for i in solution.xd
#     println(f, i)
#   end
# end 

# outfile = "t.txt"
# open(outfile, "w") do f
#     println(f, 0)
#     for i in solution.td
#         println(f, i * t_burn + t_coast)
#     end
#     println(f, t_coast + t_burn + t_coast2)
# end 