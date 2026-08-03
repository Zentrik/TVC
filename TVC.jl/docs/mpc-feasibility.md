# Why the landing guidance problem "often isn't solvable"

This is a write-up of an investigation into why the SCP landing-trajectory
problem in [`src/Guidance/6dof fixed t_burn udotdot.jl`](../src/Guidance/6dof%20fixed%20t_burn%20udotdot.jl)
frequently fails to return a trajectory when it is re-solved from the current
state, as [`Examples/MPC_Simulation.jl`](../Examples/MPC_Simulation.jl) does
every 0.25 s.

Everything below was reproduced with `PTR` and the parameters in
[`src/Guidance/run.jl`](../src/Guidance/run.jl), with `ECOS` where §1 talks
about the original behaviour and `Clarabel` afterwards.
[`Examples/FeasibilitySweep.jl`](../Examples/FeasibilitySweep.jl) reruns the
experiments.

## Summary

There are two independent problems, and they have to be separated because they
call for completely different fixes.

1. **The solve does not fail because the problem is infeasible — ECOS is
   failing.** PTR relaxes the dynamics *and* both boundary conditions with
   virtual controls, so the conic subproblem it hands to the solver is
   essentially always primal feasible. Every failure observed was
   `SCP_FAILED (NUMERICAL_ERROR)` — ECOS breaking down — or an outright
   `SingularException` crash inside the toolbox, and in the `NUMERICAL_ERROR`
   cases PTR had *already* produced a usable trajectory on an earlier iteration
   and threw it away. Handing the identical problems to Clarabel instead solves
   all of them.
2. **The problem is nonetheless nearly uncontrollable in the vertical axis
   once the motor is lit**, so even when it does solve, the MPC has almost no
   authority to correct a disturbance. This is a formulation issue, not a bug.
3. **`‖T‖ ≤ 1` was a lossy relaxation, not a lossless convexification.** The
   optimiser really does throttle down — to 89.5% on one measured case — so the
   plans were not flyable by a solid motor. Pinning `‖T‖ = 1` costs a lot of
   tractability, which is bought back by letting the touchdown time float
   instead of fixing it at burnout.

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
2. **Do not expect the PTR tolerances to save you.** Stopping before the bad
   subproblem sounds attractive, but the convergence test cannot see it coming:
   on the iteration *before* the failure the relative cost improvement is
   96–99.8% and the deviation is ~6e-2, nowhere near `ε_rel = 1e-3` or
   `ε_abs = 1e-5`. Rerunning the sweep with `ε_rel = 1e-2` changes essentially
   nothing (12/21 solved either way). The knob that works is (1): take the last
   safe iterate rather than trying to stop on one.
3. **Judge the solution on `vd`/`vbc`, not on the status string.** `SCP_SOLVED`
   only means "the last subproblem solved" — it is set even when PTR ran out
   of iterations with large virtual controls, i.e. when the returned trajectory
   does not satisfy the dynamics or the boundary conditions.
   `Examples/MPC_Simulation.jl` currently accepts any `SCP_SOLVED` plan.
4. **Change conic solver — this is the one that actually fixes it.** The README
   already notes that Mosek and Gurobi sometimes work where ECOS fails.
   [Clarabel](https://github.com/oxfordcontrol/Clarabel.jl) is open source and
   is a drop-in `solver` for `PTR.Parameters`. On the sweep above it goes from
   **12/21 to 21/21** — including the `SingularException` at `h0 = 22 m`:

   | `h0` (m) | 16 | 17 | 18 | 19 | 20 | 21 | 22 | 23 | 24 | 25 | 26–36 |
   |---|---|---|---|---|---|---|---|---|---|---|---|
   | ECOS | FAIL | ok | FAIL | FAIL | FAIL | FAIL | **crash** | FAIL | FAIL | FAIL | ok |
   | Clarabel | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok | ok |

   It also converges to a genuinely better answer where ECOS gave up: on
   `h0 = 18…25 m` Clarabel reaches `J ≈ 2e-7`, a 0.5 mm/s touchdown, against
   the `~1e-3` (3 cm/s) of the iterate ECOS discarded. The cost is iteration
   count — Clarabel took 5–47 SCP iterations against ECOS's 3–30, and it is the
   awkward low-altitude cases that need the most, so budget for that if this
   has to run at 4 Hz. Swap it in with

   ```julia
   solver, options = Clarabel, Dict("verbose" => false)
   ```

   in `src/Guidance/run.jl`.

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
one 2.5 s in can absorb about 30 cm.

Note what those numbers are measuring, though: they are the authority of a
throttle the vehicle **does not have**. The motor is solid. Aerodynamics are
deliberately left out of the guidance model (see below), so `Aero = true` in the
simulation is an unmodelled disturbance the guidance has to absorb — at 12 m/s
the axial force is ~0.25 N, ~2.4% of weight, worth ~0.8 m/s and >1 m of altitude
over a burn. That is larger than the entire throttle budget from about 2 s into
the burn, and larger still than the real budget without a throttle.

This is a divergence from the formulation in the paper
([`Paper/Guidance.tex`](../../Paper/Guidance.tex)), which uses a unit thrust
direction (correct for a solid motor) and a **free** touchdown time `b` with
`b ≥ burn time`. Fixing touchdown at burnout and adding a throttle instead is
what produced the rigidity — see §3.

### Suggested changes

* **Make the touchdown time free** (done, see §3). A second parameter `t_land`
  scales the dynamics exactly as the coast time does, bounded *below* by
  `BurnTime` so the rocket falls ballistically to the ground rather than
  arriving under thrust.
* **Fall back gracefully.** As `t0 → BurnTime` the horizon and the authority
  both go to zero and re-solving is pointless. Stop re-planning below some
  remaining-burn threshold and fly the last good plan.

Aerodynamics stay out of the guidance model on purpose — the formulation in
`Utils/Aerodynamics.jl` is not trusted, and in the simulation it is mostly there
to inject a disturbance. The consequence is just that the drag figure above is
part of the error budget rather than something the planner can anticipate.

## 3. `‖T‖ ≤ 1` is not a lossless relaxation of `‖T‖ = 1` here

The motor is solid, so the thrust magnitude is not a control — the vehicle flies
`‖T‖ = 1` whatever the plan says. The problem as written only asked for
`‖T‖ ≤ 1`.

Relaxing a nonconvex thrust bound like that is a rigorous and well known
technique — *lossless convexification*, from the GFOLD line of work (Açıkmeşe &
Ploen 2007; Blackmore, Açıkmeşe & Scharf 2010). The construction there is

```
minimise  ∫ Γ dt
s.t.      ‖T‖ ≤ Γ,   ρ_min ≤ Γ ≤ ρ_max
```

with a **slack variable Γ** that replaces `‖T‖` everywhere it appears in the cost
and in the mass dynamics. The theorem is that the optimum satisfies `‖T‖ = Γ`
pointwise, so the relaxed convex problem solves the original nonconvex one
exactly. It is also, as expected, much more tractable than the nonconvex
problem — that part is real here too, see the numbers below.

The hypothesis doing the work is that **Γ, not `‖T‖`, is what the cost sees**.
Minimising fuel drives Γ down until it meets `‖T‖`, and only then does
`ρ_min ≤ Γ` bite. This problem has none of that structure:

* there is no slack variable and no lower bound — just `‖T‖ ≤ 1`;
* the cost is `‖v_N − v_N*‖²`, which does not involve `‖T‖` at all, so nothing
  pushes the solution onto the boundary of the cone;
* the results are proved for 3-DoF translational dynamics with the thrust vector
  as the direct control. Here `T` is a *state*, rate limited by `‖Ṫ‖ ≤ 5 °/s` and
  `‖T̈‖ ≤ 10 °/s²`, and it also drives the attitude dynamics.

So it is a plain relaxation, and whether it happens to be tight is an empirical
question. It is easy to check: if the relaxation were lossless the two problems
would have the same optimal cost, and the relaxed solution would come out with
`‖T‖ = 1` anyway.

Taking `h0 = 18 m` (Clarabel, touchdown fixed at burnout in both cases):

| | `sol.cost` | `min‖T‖` | touchdown speed |
|---|---|---|---|
| `‖T‖ ≤ 1` (relaxed) | `6.98e-06` | **0.8952** | 0.0004 m/s |
| `‖T‖ = 1` (real vehicle) | `1.77e-03` | 1.0000 | 0.0267 m/s |

The relaxed optimum throttles down to **89.5%** and buys a two-orders-of-magnitude
better cost with it — 0.4 mm/s of touchdown speed against 2.7 cm/s. If the
relaxation were lossless those two rows would agree. They do not, so `‖T‖ ≤ 1` is
solving a strictly easier problem than the one the vehicle can fly, and the extra
freedom is exactly the throttle a solid motor does not have.

The relaxation *is* much more tractable, as expected — that part of the
intuition is right, and it is worth being explicit about the size of the effect:

| ignition altitude sweep (21 cases) | ECOS | Clarabel |
|---|---|---|
| `‖T‖ ≤ 1`, touchdown at burnout | 12/21 | 21/21, 5–47 SCP iterations |
| `‖T‖ = 1`, touchdown at burnout | 6/21 | solves, but hits the 50 iteration cap |

So the relaxation was buying real tractability. It just was not free.

`Examples/FeasibilitySweep.jl` prints `min‖T‖` for exactly this reason — if the
relaxed solution comes back with `min‖T‖ < 1` it is planning a throttle the
vehicle does not have, and the simulation's `normalize(...)` will quietly throw
that part of the plan away.

### What changed

`RocketParameters` gained two switches, both defaulting to the physical vehicle:

* `Throttleable = false` adds `‖T‖ ≥ 1` as a nonconvex path constraint,
  linearised about the reference by the SCP algorithm (`s(x) = 1 − T·T ≤ 0`) and
  relaxed by a penalised virtual control, so it can never make a subproblem
  infeasible. Together with the existing cone this pins `‖T‖ = 1`. Set it to
  `true` to get the old relaxation back.
* `FixedLandingTime = false` promotes the touchdown time to a decision variable
  `p[veh.id_tland] ∈ [BurnTime, BurnTime + MaxBallisticTime]`, and the
  trajectory is scaled to that horizon rather than to burnout. Past `BurnTime`
  the thrust table returns 0 and the mass and CG tables are flat, so the same
  dynamics carry straight on as an unpowered fall — no extra modelling needed.

  Touchdown is **not** allowed before burnout. Thrust to weight is ~1.35 for
  most of the burn, so a rocket that reaches the ground while still thrusting
  bounces and flies again; the paper counts that as a failed landing, and
  `height ≥ 0` at every node is what holds the trajectory up until the motor is
  spent.

The second is what pays for the first. Requiring touchdown *exactly* at burnout
was only tractable because the throttle was there to absorb the terminal
altitude constraint; with the throttle gone the constraint has nothing to work
with, which is what §2 is about. Allowing a ballistic tail turns "be exactly at
the ground the instant the motor cuts out" into "be above the ground at cutout
and fall the rest of the way" — the condition the paper's `b ≥ burn time`
expresses, physically meaningful, and unlike the throttle something the vehicle
can actually deliver.

On the nominal ignition state the two switches together give

```
|T|=1, free t_land   SCP_SOLVED   t_coast=1.327  t_land=3.450
                                  min‖T‖=1.0000  r_end=[32.6, -13.44, 0.0]  |v_end|=0.228 m/s
```

against 0.52 m/s for the old relaxed-throttle, fixed-touchdown problem — a
better landing, and one the vehicle can actually fly. Note that the optimiser
picked `t_land = BurnTime` here: on the nominal trajectory the freedom is not
needed, it is there for when a disturbance means the old problem would have had
no answer at all.

Restarting part way through the burn — the case the MPC actually depends on,
and the one that used to throw `SingularException` on all 20 attempts — now
solves every time:

| restart | `t_land` | `min‖T‖` | `h_end` | touchdown |
|---|---|---|---|---|
| `t0 = 1.5 s`, on the nominal | 3.450 | 1.0000 | 0.000 | 0.52 m/s |
| `t0 = 1.5 s`, 0.5 m low | 3.450 | 1.0000 | 0.000 | 1.31 m/s |
| `t0 = 1.5 s`, 0.5 m high | 3.450 | 1.0000 | 0.000 | 1.10 m/s |
| `t0 = 1.5 s`, 2 m low | **3.784** | 1.0000 | 0.000 | 2.72 m/s |
| `t0 = 1.5 s`, 2 m high | 3.450 | 1.0000 | 0.000 | 2.82 m/s |
| `t0 = 2.5 s`, on the nominal | 3.450 | 1.0000 | 0.000 | 0.52 m/s |
| `t0 = 2.5 s`, 0.5 m low | 3.450 | 1.0000 | 0.000 | 0.52 m/s |
| `t0 = 2.5 s`, 0.5 m high | **3.720** | 1.0000 | 0.000 | 3.17 m/s |
| `t0 = 2.5 s`, 2 m low | 3.450 | 1.0000 | 0.000 | 0.52 m/s |
| `t0 = 2.5 s`, 2 m high | **4.037** | 1.0000 | 0.000 | 6.28 m/s |

10/10 solved, every one touching down at zero altitude with `‖T‖ = 1` — no
phantom throttle anywhere. Most cases still land exactly at burnout, which is
the right answer when the state allows it. The bolded rows are the ones that
use the ballistic tail: 0.27–0.59 s of unpowered fall, which is what makes a
state that cannot reach the ground by burnout solvable at all instead of
leaving the terminal altitude constraint with no answer.

**Not yet measured**: a full ignition-altitude sweep with both switches on.
`Examples/FeasibilitySweep.jl` runs it.

## 4. Closed loop results

[`Examples/MPCSweep.jl`](../Examples/MPCSweep.jl) flies the whole thing: TVC's
own `f!` as the plant with `Aero = true` (the disturbance the guidance never
sees), thrust direction taken from the current plan and applied at full motor
thrust, and `solveProblem` re-run every 0.25 s from the measured state. Six
release states, ordered by how hard the landing was:

| release state | touchdown `‖v‖` | `v_z` | `v_xy` | tilt | `‖ω‖` | contact at motor time | solves |
|---|---|---|---|---|---|---|---|
| 5 m lower | 1.74 | −1.53 | 0.84 | 0.8° | 10.04 | 3.39 (**before** burnout) | 15, 0 rejected |
| nominal | 2.46 | −2.39 | 0.60 | 0.5° | 0.15 | 3.45 (at burnout) | 16, 0 rejected |
| already descending 3 m/s | 2.92 | −2.85 | 0.65 | 0.5° | 0.03 | 3.45 (at burnout) | 15, 0 rejected |
| released tilted 5° | 4.64 | −4.10 | 2.18 | 11.9° | 1.43 | 3.41 (**before** burnout) | 16, 0 rejected |
| faster horizontal (6, −5) | 9.03 | −7.99 | 4.21 | 16.0° | 4.62 | 2.86 (**before** burnout) | 14, 0 rejected |
| 5 m higher | 12.71 | −12.12 | 3.83 | 33.5° | 0.42 | 4.27 (0.82 s ballistic) | 15, 0 rejected |

**The solver is no longer the problem.** 91 guidance solves across six flights,
zero rejected — no `SCP_FAILED`, no `SingularException`, nothing thrown. That
was the original complaint and it is gone.

**What is left is control authority, and it is the vertical channel.** Look at
`5 m higher`: it burns out roughly 6.6 m up still moving down at ~4 m/s, then
falls the remaining distance, arriving at 12.1 m/s. That is not the optimiser
giving up — with `‖T‖ = 1` and a fixed impulse there is no trajectory that puts
the vehicle at the ground with low speed from that release point. The ballistic
tail made the problem *solvable*; nothing can make it *soft*. The two cases that
land well (nominal and already-descending, 2.5–2.9 m/s and 0.5° of tilt) are the
ones whose release state happens to sit near the reachable set §2 describes.

**Three of six contact the ground before burnout**, which is the bounce case:
`height ≥ 0` holds in the *plan*, but the vehicle, flying with aerodynamics the
guidance does not model and only re-planned every 0.25 s, arrives early anyway.
`faster horizontal` is 0.59 s early at 9 m/s and 16° of tilt — a crash, not a
landing.

**One thing to look at that is not authority.** Touchdown `‖ω‖` is 10.0 rad/s on
`5 m lower` and 4.6 on `faster horizontal`, on flights that are otherwise
upright — the guidance constrains `ω = 0` at the terminal node, so something is
diverging near the end. The prime suspect is roll: `I_zz = 2.48e-4` against
`I_xx = I_yy = 8.27e-2`, and the roll torque limit of 0.1 N·m buys ~400 rad/s²,
so a roll command sampled from a plan every 0.25 s has enormous rate error for
very little timing mismatch. This has not been confirmed — the sweep reports
`‖ω‖` and not its components. Worth instrumenting before reading anything else
into it.

## 5. Bugs found

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
