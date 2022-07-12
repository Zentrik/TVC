using TVC, DifferentialEquations, Plots, LinearAlgebra

#   Specify Parameters
#   ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡

veh = RocketParameters()
atmos = Atmosphere()
traj = RocketTrajectoryParameters();

mdl = RocketProblem(veh, atmos, traj)

servoΔt = 0.02 # Servo step rate
x₀ = [[0; 0; 15.]; [0; 0; 0]; [1; 0; 0; 0]; [0; 0; 0]][[veh.id_r; veh.id_v; veh.id_quat; veh.id_ω]];

#   Specify Controller
#   ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡

#   PID Controller
#   ================

function control(x, p, t)
    veh = p.veh
    tₘ = motorTime(t, p.MotorIgnitionTime)

    desired_torque = -.5 * x[veh.id_quat[2:4]] - .5 * x[veh.id_ω] # PD controller, not well tuned.

    return zeros(3)
end

#   Continuous Actuator Model
#   ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡

function Actuator(x, p, t, desired_torque) # Model of TVC
    veh = p.veh
    tₘ = motorTime(t, p.MotorIgnitionTime)

    Thrustxy = [0 1; -1 0] * desired_torque[1:2] / veh.MomentArm(tₘ)[3];
    # r \times Thrust = Mb, if r = [0; 0; z] then Thrust = (Mb(2) / z, - Mb(1) / z, 0), then we add in z component of Thrust
    # (Mb(2) / z, - Mb(1) / z) = [0 1; -1 0] * Mb / z
    # See StateSpaceController.m for context

    if veh.Thrust(tₘ) ≈ 0
        Thrust = zeros(3)
    else
        if norm(Thrustxy) < veh.Thrust(tₘ) * sind(5)
            Thrustz = sqrt(veh.Thrust(tₘ)^2 - norm(Thrustxy)^2);
            Thrust = [Thrustxy; Thrustz]
        else 
            Thrustxy = normalize(Thrustxy) * veh.Thrust(tₘ) * sind(5)
            Thrustz = veh.Thrust(tₘ) * cosd(5)
            Thrust = [Thrustxy; Thrustz]
        end
    end

    roll = clamp(desired_torque[3], -.1, .1)
    torque = veh.MomentArm(tₘ) × Thrust + [0; 0; roll]
    
    return (force=Thrust, torque=torque)
end

#   Solve ODE
#   ≡≡≡≡≡≡≡≡≡≡≡

tspan = (0, veh.BurnTime)
p = (veh=veh, atmos=atmos, Aero=false, wind=zeros(3), MotorIgnitionTime=0., Control=(x, p, t) -> Actuator(x, p, t, control(x, p, t)))

prob = ODEProblem(f!, x₀, tspan, p)

condition(x,t,integrator) = x[3] # when height is zero halt integration
cb = ContinuousCallback(condition, nothing, terminate!) # when going upwards do nothing

sol = DifferentialEquations.solve(prob, reltol=1e-8, abstol=1e-8, callback=cb)

plot(sol, vars=veh.id_r)
# plot(sol, vars=veh.id_quat)
# plot(sol, vars=veh.id_ω)

# deviationAngle = map(q -> acosd(TVC.Utils.rotate(q, [0; 0; 1]) ⋅ [0; 0; 1]), eachcol(sol[veh.id_quat, :]))

# plot(sol.t, deviationAngle)