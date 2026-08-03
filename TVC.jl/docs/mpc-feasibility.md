# Why the landing guidance problem "often isn't solvable"

This is a write-up of an investigation into why the SCP landing-trajectory
problem in [`src/Guidance/6dof fixed t_burn udotdot.jl`](../src/Guidance/6dof%20fixed%20t_burn%20udotdot.jl)
frequently fails to return a trajectory when it is re-solved from the current
state, as [`Examples/MPC_Simulation.jl`](../Examples/MPC_Simulation.jl) does
every 0.25 s.

Everything below was reproduced with `ECOS`, `PTR`, and the parameters in
[`src/Guidance/run.jl`](../src/Guidance/run.jl).
[`Examples/FeasibilitySweep.jl`](../Examples/FeasibilitySweep.jl) reruns the
experiments.

## Summary

There are two independent problems, and they have to be separated because they
call for completely different fixes.

1. **The solve does not fail because the problem is infeasible.** PTR relaxes
   the dynamics *and* both boundary conditions with virtual controls, so the
   conic subproblem it hands to ECOS is essentially always primal feasible.
   Every failure observed was `SCP_FAILED (NUMERICAL_ERROR)` — ECOS breaking
   down — or an outright `SingularException` crash inside the toolbox. In the
   `NUMERICAL_ERROR` cases PTR had *already* produced a perfectly usable
   trajectory on an earlier iteration and threw it away.
2. **The problem is nonetheless nearly uncontrollable in the vertical axis
   once the motor is lit**, so even when it does solve, the MPC has almost no
   authority to correct a disturbance. This is a formulation issue, not a bug.

Plus a handful of concrete bugs (listed at the end) that make (1) much worse.

## 1. The failures are numerical, not infeasibility

### The evidence

Sweeping only the ignition altitude and leaving everything else at the
defaults (`r0 = [20, -4, h]`, `v0 = [4, -3, 0]`, everything else zero):

| `h0` (m) | result | iterations | best iterate PTR had |
|---|---|---|---|
| 16 | `SCP_FAILED (NUMERICAL_ERROR)` | 3 | #2, `J = 7.07e-3`, `max\|vd\| = 4.1e-3` |
| 17 | `SCP_SOLVED` | 14 | #14, `J = 4.13e-6` |
| 18 | `SCP_FAILED (NUMERICAL_ERROR)` | 3 | #2, `J = 1.67e-3`, `max\|vd\| = 2.8e-3` |
| 19 | `SCP_FAILED (NUMERICAL_ERROR)` | 3 | #2, `J = 1.14e-3` |
| 20 | `SCP_FAILED (NUMERICAL_ERROR)` | 3 | #2, `J = 8.22e-4` |
| 21 | `SCP_FAILED (NUMERICAL_ERROR)` | 3 | #2, `J = 6.57e-4` |
| 22 | **`SingularException(3)` thrown** | – | – |
| 23 | `SCP_FAILED (NUMERICAL_ERROR)` | 3 | #2, `J = 4.08e-4` |
| 24 | `SCP_FAILED (NUMERICAL_ERROR)` | 3 | #2, `J = 3.31e-4` |
| 25 | `SCP_FAILED (NUMERICAL_ERROR)` | 3 | #2, `J = 2.82e-4` |
| 26 | `SCP_SOLVED` | 5 | #5, `J = 5.9e-12` |
| 27 | `SCP_SOLVED` | 7 | #7, `J = 3.7e-11` |
| 28 | `SCP_SOLVED` | 30 | #30, `J = 1.11e-2` |
| 29–33 | `SCP_SOLVED` | 5–9 | `J = 0.099 … 1.31` |

Two things stand out. The failures are *not* the hard problems — the reported
cost of the discarded iterate on the failing cases (`1e-4` to `7e-3`, i.e. a
1–8 cm/s touchdown) is far better than that of several cases that succeed
(`h0 = 30 m` converges happily to `J = 0.274`, a 0.52 m/s touchdown). And the
iteration count is erratic: 3, 14, 3, 3, …, 5, 7, 30, 9. Neither pattern is
what an infeasible or a marginally feasible problem looks like; both are what
a conic solver falling over looks like.

(This is with the three bugs at the bottom of this page fixed. Before the
fixes the same sweep failed at 18, 20, 22, 24, 28 and succeeded at 26, 30, 32,
34, which is less monotone but has the same character. The bugs are worth
fixing, but they are not what makes the solve fail.)

### What actually happens

Look at the PTR iteration log for a failing case (`h0 = 18 m`):

```
k | status   | vd    | vs    | vbc   | J        | ΔJ %  | Δx    | ... | dyn | ηx   | ηu   | ηp
1 | ALMOST_O | 1e-01 | 0e+00 | 2e-10 | 9.09e-01 |       | 1e-01 |     | T   | 0.12 | 0.40 | 0.12
2 | ALMOST_O | 3e-03 | 0e+00 | 7e-09 | 2.11e-02 | 97.67 | 5e-02 |     | T   | 0.07 | 0.61 | 0.02
unsafe solution (NUMERICAL_ERROR), exiting
```

Iteration 2 is already a good answer: the dynamics virtual control is 3e-3,
the boundary-condition virtual control is 7e-9 (i.e. the terminal conditions
*are* met), and it is dynamically feasible. PTR then takes a third step, ECOS
returns `NUMERICAL_ERROR`, and `SCPToolbox` throws the whole thing away:

```julia
# SCPToolbox/src/solvers/scp.jl, SCPSolution(history)
if unsafe_solution(last_sol)      # <- only ever looks at the LAST subproblem
    status = SCP_FAILED ...
    xd = RealMatrix(undef, ...)   # <- previous iterations discarded
    cost = Inf
```

Note also that **ECOS returns `ALMOST_OPTIMAL`, not `OPTIMAL`, on nearly every
iteration of every solve**, including the ones that succeed. The subproblem is
permanently sitting on ECOS's accuracy limit; whether a particular initial
condition "works" comes down to whether PTR's stopping criterion happens to
trigger before ECOS tips over. That is why the behaviour looks random with
respect to the initial state.

Why the subproblem is so badly conditioned:

* The trust region is soft (`q_tr = Inf`, penalty `wtr * η`). As PTR converges
  the radii `ηx, ηu, ηp` go to 0, so the constraints `|x̂ - x̂_ref| ≤ η`
  degenerate into equalities and the subproblem loses strict feasibility.
  The failing iteration is always the one after the iterate has essentially
  converged.
* Several second-order cones are active simultaneously at the solution
  (`‖T‖ ≤ 1` and `‖q‖ ≤ 1` are both tight along most of the trajectory, and
  the gimbal cone is tight wherever the vehicle is steering).
* The roll axis is scaled ~300x differently from the pitch/yaw axes
  (`I_zz = 2.48e-4` against `I_xx = I_yy = 8.27e-2`), so `B[ω, roll]` is ~4000
  while the rest of the input matrix is O(1)–O(10).
* The terminal cost `dot(x_v - vN, x_v - vN)` is a *quadratic* objective, and
  ECOS is a pure conic solver, so MOI bridges it into an extra rotated
  second-order cone. `SCPToolbox`'s own `QuadraticCost` docstring says "for
  robustness (in JuMP) it has been observed that it is best to reformulate the
  problem (via epigraph form) such that this function is affine". Minimising
  `‖v_N - vN‖` against an explicit epigraph variable, rather than the squared
  norm, is the same minimiser and one fewer bridge.

### What to do about it

In rough order of value for effort:

1. **Never discard a good iterate.** For MPC the right output is "the last
   iterate whose subproblem solved", not "nothing". The information is already
   in the history that `PTR.solve` returns:

   ```julia
   sol, hist = PTR.solve(ptr_pbm)
   if startswith(sol.status, "SCP_FAILED")
       subs = [s.sol for s in hist.subproblems]
       k = findlast(s -> s.status == MOI.OPTIMAL || s.status == MOI.ALMOST_OPTIMAL, subs)
       # subs[k].xd / .ud / .p is a usable trajectory; check subs[k].feas and
       # max(abs, subs[k].vd) before accepting it
   end
   ```

   Every failing case inspected had one: e.g. `h0 = 16 m` reports
   `SCP_FAILED (NUMERICAL_ERROR)` but iteration 2 has `J = 7.07e-3` and
   `max|vd| = 4.1e-3`, which is a perfectly good landing trajectory.
2. **Stop iterating earlier.** The failing solves fail on the iteration *after*
   they have effectively converged, so a looser `ε_rel` (`1e-2` rather than
   `1e-3`) should avoid a good fraction of them, and for MPC there is no value
   in the last 1% of cost. Worth measuring on your machine — the failure point
   depends on the exact ECOS build.
3. **Judge the solution on `vd`/`vbc`, not on the status string.** `SCP_SOLVED`
   only means "the last subproblem solved" — it is set even when PTR ran out
   of iterations with large virtual controls, i.e. when the returned trajectory
   does not satisfy the dynamics or the boundary conditions.
   `Examples/MPC_Simulation.jl` currently accepts any `SCP_SOLVED` plan.
4. **Try a different conic solver.** The README already notes that Mosek and
   Gurobi sometimes work where ECOS fails. Clarabel is worth a look too — it
   is open source, handles the SOC/exponential cones this problem uses, and is
   generally better behaved than ECOS on badly scaled problems.

## 2. The vehicle has almost no vertical authority once the motor is lit

This is why the MPC cannot do its job even when the solver behaves.

The problem as coded requires touchdown to happen **exactly at motor burnout**:
the last node is at motor time `BurnTime`, and the terminal boundary condition
is `r_z = 0` there. The knobs available to hit that one scalar condition are:

* the coast time before ignition (`p[1]`) — **hard-fixed to 0 as soon as
  `MotorFired` is true**;
* the "throttle" `‖T‖ ≤ 1` — but `‖Ṫ‖ ≤ 5 deg/s` limits how fast `‖T‖` can
  move, and it starts at 1;
* tilting the whole vehicle, which costs `cos θ` of vertical thrust but is
  coupled to the horizontal channel and to the terminal attitude constraint.

Integrating the vertical channel with the real F15 curve (49.6 N·s over 3.45 s,
mass 1.093 → 1.033 kg, so T/W ≈ 1.4 for most of the burn: at full throttle the
vehicle is always decelerating, never coasting down) gives the total throttle
authority over a whole burn:

| | full throttle | minimum throttle | span |
|---|---|---|---|
| Δh over the burn | +24.09 m | +17.43 m | **6.66 m** |
| Δv over the burn | +12.82 m/s | +6.96 m/s | **5.85 m/s** |

and, more importantly, what is left of it partway through the burn:

| burn remaining | 3.45 s | 2.95 s | 2.45 s | 1.95 s | 1.45 s | 0.95 s | 0.45 s |
|---|---|---|---|---|---|---|---|
| altitude window | 6.66 m | 3.89 m | 2.08 m | 0.95 m | 0.34 m | 0.07 m | 0.004 m |
| velocity window | 5.85 m/s | 4.09 m/s | 2.68 m/s | 1.57 m/s | 0.77 m/s | 0.26 m/s | 0.03 m/s |

So a re-solve 2 s into the burn can absorb under a metre of altitude error, and
one 2.5 s in can absorb about 30 cm. Two things guarantee that much error:

* **The guidance model has no aerodynamics at all**, while the simulation runs
  with `Aero = true`. At 12 m/s the axial force is ~0.25 N, ~2.4% of weight,
  which is ~0.8 m/s and >1 m of altitude over a burn — larger than the entire
  correction budget from about 2 s in. The airframe has no fins, so the normal
  force also produces a pitching moment the guidance never sees.
* **The controller does not fly the plan.** `Examples/MPC_Simulation.jl` does

  ```julia
  desired_tvc = normalize(sample(sol.xc, time)[veh.id_T])   # magnitude thrown away
  ...
  Thrust = desired_tvc * veh.Thrust(tₘ)                     # always 100%
  ```

  so whenever the guidance plans `‖T‖ < 1` the vehicle flies full thrust
  instead. That is *the* vertical control channel being discarded, and the
  error it injects compounds every 0.25 s.

This is also a divergence from the formulation in the paper
([`Paper/Guidance.tex`](../../Paper/Guidance.tex)), which uses a unit thrust
direction (no throttle at all, correct for a solid motor) and a **free**
touchdown time `b` with `b ≥ burn time`. Fixing touchdown at burnout and adding
a throttle the vehicle does not have is what produced the rigidity.

### Suggested changes

* **Make the touchdown time free.** Add a second parameter `t_land` and scale
  the dynamics by it, exactly as the coast time is handled now, with
  `t_land ≤ BurnTime`. This restores the free variable the paper's formulation
  has and gives the mid-burn re-solve something to move.
* **Either honour the throttle or remove it.** If the motor cannot throttle,
  constrain `‖T‖ = 1` (linearised about the reference: `T_ref·T = 1`) and let
  vehicle tilt be the vertical knob. If a partial throttle is genuinely
  available, apply the planned magnitude in the simulation instead of
  normalising it.
* **Put the drag term in the guidance dynamics**, at least the axial component.
  It is a smooth function of the state and cheap to linearise, and it is
  currently larger than the correction authority it is competing with.
* **Fall back gracefully.** As `t0 → BurnTime` the horizon and the authority
  both go to zero and re-solving is pointless. Stop re-planning below some
  remaining-burn threshold and fly the last good plan.

## 3. Bugs found

### `slerp_quat` ignores its interpolation parameter (fixed)

`src/Utils/Quaternions.jl`:

```julia
axis, angle = quatLogAxisAngle(Δq)
Δq_t = [cos(angle); sin(angle) * axis]     # `frac` is never used
```

`slerp_quat(q0, q1, frac)` returned `normalize(q1)` for *every* value of
`frac`. It is used in exactly one place — the initial guess:

```julia
x[veh.id_quat, k] = slerp_quat(traj.q0, [traj.q0[1]; traj.qN[2:3]; traj.q0[4]], mix)
```

so the attitude initial guess was a single constant quaternion, equal to the
*target* attitude, at every node, and **the vehicle's current attitude never
entered the guess at all**. That is the worst possible reference for the first
linearisation, and it is worst precisely in the MPC case, where the vehicle is
tilted and rotating. It is also the direct cause of the note at the top of
`6dof fixed t_burn udotdot.jl` about the reference trajectory having non-unit
quaternions: `[q0[1]; qN[2:3]; q0[4]]` has norm < 1 whenever the rocket is
tilted, and it was never normalised.

### `wexp` is not differentiable at zero, and silently returned a zero Jacobian (fixed)

`src/Utils/Quaternions.jl` used `theta = norm(w)` with a `theta < eps()` early
return of a constant. `norm` is not differentiable at `w = 0` (it produces a
`NaN` partial), and the early return means `ForwardDiff` sees a *constant* and
reports a derivative of zero. The initial-condition Jacobian with respect to
the coast time does exactly this:

```julia
J[veh.id_quat, veh.id_tcoast] = ForwardDiff.derivative(t -> quatL(traj.q0) * wexp(traj.ω0 * t), p[veh.id_tcoast])
```

and the reference value of `p` is 0 both for the default initial guess and
every time `MotorFired` is true. So the first PTR iteration always believed
that coasting does not rotate the vehicle. With `ω0 ≈ 0.5 rad/s` — typical in
the saved states in `Examples/FeasibilityTests.jl` — a 0.4 s coast rotates the
airframe by more than 10°, and the linearisation said zero.

`wexp` is now written in terms of `θ² = w·w` with a series expansion near zero,
which is smooth everywhere and agrees with the old expression to machine
precision away from zero.

### State 13 scaling typo (fixed)

```julia
advise!(pbm, :state, 11, (-10.0, 10.0))
advise!(pbm, :state, 12, (-10.0, 10.0))
advise!(pbm, :state, 13, (-10.0, 00.0))   # ω_z is not sign definite
```

`compute_scaling` turns this into `ω_z = 10 ω̂_z - 10`, i.e. `ω_z = 0` maps to
`ω̂_z = 1` and the scale factor is half that of the other two axes.

### `SingularException` crash instead of a solver failure

Every mid-burn re-solve tried (20/20 restarts off the nominal trajectory, at
`t0 ∈ {0.5, 1.5, 2.5, 3.0}` s with altitude errors between −2 m and +2 m), and
all five states saved in `Examples/FeasibilityTests.jl`, threw
`SingularException(1)` out of `PTR.solve`, from

```
SCPToolbox/src/solvers/discretization.jl:267   iPhi = Phi \ I(nx)
```

reached via `SubproblemSolution(spbm)` → `discretize!`. Instrumenting the
toolbox shows what is actually going on:

```
>>> iteration 1: termination_status = ALMOST_OPTIMAL, primal_status = NEARLY_FEASIBLE_POINT,
                 result_count = 1, maxabs(x) = NaN, maxabs(u) = NaN, p = [NaN]
!!! SINGULAR Phi at t=0.00056 k=1
  x = [Inf, Inf, Inf, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, Inf, Inf, Inf, NaN, NaN, NaN]
  u = [Inf, NaN, Inf, NaN]
```

**ECOS reports `ALMOST_OPTIMAL` / `NEARLY_FEASIBLE_POINT` and hands back a
solution vector that is entirely `NaN`.** `SCPToolbox`'s `unsafe_solution`
accepts `ALMOST_OPTIMAL` as a good solve, so the `NaN`s are taken as a
trajectory, `discretize!` propagates a `NaN` state transition matrix, and the
`lu` inside `Phi \ I(nx)` reports it as singular. Two things follow:

* `solveProblem` can throw rather than return a status, and
  `Examples/MPC_Simulation.jl` has no `try`/`catch` around it, so this aborts
  the whole simulation instead of skipping one MPC tick.
* Even without the crash, an `ALMOST_OPTIMAL` result is not necessarily
  usable. Check `all(isfinite, sol.xd)` before accepting a plan, and consider
  treating `ALMOST_OPTIMAL` as suspect — note from the logs above that ECOS
  returns it on essentially *every* iteration of *every* solve here, including
  the ones that converge.

### No guard for `t0 ≥ BurnTime`

`mpc!` fires on a `PeriodicCallback` over
`tspan = (t0, t0 + coast + BurnTime + 10)`, with no check that the motor still
has burn time left. Once `t0 > BurnTime`, `BurnTime - t0 < 0`, so the dynamics
are scaled by a negative number (time runs backwards), `N` clamps to its floor
of 5, and `veh.Thrust` has already run off the end of its table and returns 0,
leaving the thrust states completely uncontrollable.

### Minor

* `mpc!` reads `p.solution[]` and `p.t0[]` from the global `p` rather than
  `integrator.p`. They happen to be the same object, so it works today.
* The coast-phase initial condition assumes `ω` is constant during free-fall.
  There is already a `println` warning for this; with `I_xx = I_yy` it is only
  wrong to the extent that aerodynamic moments act, which for a finless
  airframe at ~10 m/s is not nothing.
