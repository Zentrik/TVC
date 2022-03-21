using XLSX, Interpolations

ThrustTable = XLSX.readdata("$(@__DIR__)/../data/F15_Thrust.xlsx", "Sheet1!A1:B27"); # $(@__DIR__) gives directory of Rocket_Acceleration.jl
MassTable = XLSX.readdata("$(@__DIR__)/../data/mass,cg,I.xlsx", "mass,cg,I!A2:B3608");
CGTable = [XLSX.readdata("$(@__DIR__)/../data/mass,cg,I.xlsx", "mass,cg,I!A2:A3608") XLSX.readdata("$(@__DIR__)/../data/mass,cg,I.xlsx", "mass,cg,I!E2:E3608")]

ThrustInterp = LinearInterpolation(ThrustTable[:, 1], ThrustTable[:, 2], extrapolation_bc=0)
MassInterp = LinearInterpolation(MassTable[:, 1], MassTable[:, 2], extrapolation_bc=Flat())
CG = LinearInterpolation(CGTable[:, 1], CGTable[:, 2], extrapolation_bc=Flat())

Thrust(t) = ThrustInterp(t) 
Mass(t) = MassInterp(t)
Acceleration(t) = Thrust(t) / Mass(t)

function ThrottleLevel(u_norm, norm_max)
    len = length(u_norm)
    Throttle = zeros(len)

    for k = 1:len
        if u_norm[k] <= 1e-4 && norm_max[k] <= 1e-4 # prevents divide by 0 and accounts for floating point error?
            Throttle[k] = 1
        else
            Throttle[k] = u_norm[k] / norm_max[k]
        end
    end
    return Throttle
end