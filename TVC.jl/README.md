# TVC


## Notes
SCP Problem requires 2 iterations to stop, as we need to see change in solution between two solves. (1st iteration is not the initial guess but the iteration from it).

SCS/ ProxSDP/ COSMO dont seem to work well
To try: COPT/ SeDuMi/ SDPT3
Mosek and Gurobi work ok, they seem to work when ECOS may fail but don't converge/ slow to converge in most cases especially when ECOS works fine.
ECOS works decently

## To Do
- Change `I(n)` to `I`.
- More accurate actuator model, account for torque exerted by gimbal moving.
- @unpack?
- Switch to Rotations.jl
- ModellingToolkit.jl? 
<br>
<br>
- Add realistic wind.
<br>
<br>
- Look at style guide and fix capital letters/_ in function/variable struct names.
<br>
<br>
- Deal with when lcvx won't hold, i.e. when intersection of optimal trajectories and interior of reachable set is non empty? 
Basically need to add targeting of minimum landing position error if thrust constraint is violated because it can easily land with 0 velocity.
- Allow for rocket to land after motor finishes burning.
- Allow rocket to land earlier when close to surface/ motor about to finish otherwise it might be infeasible?
- Penalise horizontal velocity more when landing than vertical to prevent tipping over?
<br>
<br>
<br>
<br>
- Add MEKF (https://github.com/ChiyenLee/EKF.jl)

### Performance Improvements
- Deal with functions in structs better, e.g. use function wrappers or parameterise their type.
- `@deftype Union{SVector{3, R}, Vector{R}}`, is this fine or should it be parameterised?
- Use StaticArrays more?