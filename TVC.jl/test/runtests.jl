using TVC
using Printf
using SCPToolbox
using Test

@testset "Guidance.jl" begin
    # Write your tests here.'
    @test solve().status == @sprintf("%s", SCP_SOLVED)

    veh = RocketParameters()
    atmos = Atmosphere()
    traj = RocketTrajectoryParameters(MotorFired=true);

    mdl = RocketProblem(veh, atmos, traj)
    @test solve(mdl).status == @sprintf("%s", SCP_SOLVED)
end