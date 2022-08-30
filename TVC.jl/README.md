# TVC


## Notes
SCP Problem requires 2 iterations to stop, as we need to see change in solution between two solves. (1st iteration is not the initial guess but the iteration from it).

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
<br>
<br>
- Add MEKF (https://github.com/ChiyenLee/EKF.jl)

### Performance Improvements
- Deal with functions in structs better, e.g. use function wrappers or parameterise their type.
- `@deftype Union{SVector{3, R}, Vector{R}}`, is this fine or should it be parameterised?
- Use StaticArrays more?