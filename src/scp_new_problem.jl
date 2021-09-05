module scp_new_problem
include("./3dof_testing.jl")
end # module

using .scp_new_problem
using Utils, LinearAlgebra, Plots

solution = scp_new_problem.solve(); # Remember J is augmented cost function
# print(sample(solution.xc, 1)[4:6]' * sample(solution.xc, 1)[4:6])

ctres = 1000
N = size(solution.xd, 2)
xct = hcat([sample(solution.xc, t) for t in LinRange(0, 1, ctres)]...)
vct = hcat([sample(solution.uc, t) for t in LinRange(0, 1, ctres)]...)

plotlyjs()

t_coast = 3.45 # solution.p[2]

Plots.plot(LinRange(0, 1, ctres) .* t_coast, xct[1:3, :]', title = "Position")

Plots.plot(LinRange(0, 1, ctres) .* t_coast, xct[4:6, :]', title = "Velocity")

Plots.plot(LinRange(0, 1, ctres) .* t_coast, vct')
Plots.plot(LinRange(0, 1, ctres) .* t_coast, norm.(eachcol(vct)), title = "Thrust Magnitude") # norm constraint is broken intersample ...

println(solution.p[1])

v_N =[0; 0; 0];
println(solution.cost^0.5)


# using Solvers
# fieldnames(Solvers.SCPSolution)
