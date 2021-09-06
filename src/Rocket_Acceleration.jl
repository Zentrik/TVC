using XLSX, Interpolations

ThrustTable = XLSX.readdata("E:/Code/TVC/Rocket/F15_Thrust.xlsx", "Sheet1!A1:B27");
MassTable = XLSX.readdata("E:/Code/TVC/Rocket/mass,cg,I.xlsx", "mass,cg,I!A2:B3608");

ThrustInterp = LinearInterpolation(ThrustTable[:, 1], ThrustTable[:, 2], extrapolation_bc=0)
MassInterp = LinearInterpolation(MassTable[:, 1], MassTable[:, 2], extrapolation_bc=Flat())
Acceleration(t) = ThrustInterp(t) ./ MassInterp(t)

function ThrottleLevel(u_norm, norm_max)
    len = length(u_norm)
    Throttle = zeros(len)

    for k = 1:len
        if u_norm[k] <= 1e-6 && norm_max[k] <= 1e-6 # prevents divide by 0 and accounts for floating point error?
            Throttle[k] = 1
        else
            Throttle[k] = u_norm[k] / norm_max[k]
        end
    end
    return Throttle
end