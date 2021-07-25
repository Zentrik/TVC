module scp_new_problem
include("./my_problem.jl")
end # module

using .scp_new_problem

solution = scp_new_problem.solve()
scp_new_problem.plot(solution)
