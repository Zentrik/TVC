module TVC

function motorTime end

include("Utils/Utils.jl")
using .Utils 
export Atmosphere, RocketParameters, f!


include("Guidance/Guidance.jl")
using .Guidance
export solveProblem, RocketProblem, RocketTrajectoryParameters, printSolution, plotThrottle, saveTrajectory

export motorTime

end