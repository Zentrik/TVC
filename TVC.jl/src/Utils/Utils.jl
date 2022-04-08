module Utils

include("Rocket_Acceleration.jl")

include("Quaternions.jl")
export quatL, quatR, skew, wexp, wexp_w, conjugate, rotate, slerp, slerp_quat, qualLog, to_matrix

include("Aerodynamics.jl")
include("Atmosphere.jl")
include("Parameters.jl")

include("Dynamics.jl")

end