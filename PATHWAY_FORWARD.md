# Pathway forward

Working plan, ordered by what blocks what. Updated 2026-08-13, at v1.12.0.

This is a live document: close items by deleting them, and record what was
learned in the place that will be read again (`PROJECT_RULES.md` for a rule,
`reference/<bench>/README.md` for a benchmark finding, a commit message for a
fix). Do not let this file become the archive.

---

## 0. State of play

Multi-fault support is real but young. The engine, the mesh generator and the
on-fault input path are all neutral to `ntotft`, and a two-fault step-over runs
end to end. Three multi-fault bugs have been found so far — accumulator
aliasing, uninitialised output slabs, and a fault plane falling between mesh
lines — and **every one of them was found by running a new case, not by the
test suite**. That is the single most important fact on this page. See item 3.

BP8 conforms to the 2026-08-12 benchmark description. The over-determined
initial condition that cost days is resolved, in favour of what we already
shipped. What remains open on BP8 is ours, not the spec's: the results are not
converged in domain size.

---

## 1. Submission — BP8-QD-GS and BP8-QD-PW at 50 m

Ship the coarse pair for a quick platform view, with the convergence caveat
stated. Steps: resample profiles onto the 81-node 10 m grid required by section
4.3 (`scripts/resampleBP8Profiles.py`), validate (`scripts/checkBP8Submission`),
package as `dliu_eqquasi-<version>-50m.zip`.

The 10 m production entry waits on item 2. There is no point spending 25 hours
of wall clock refining a mesh whose size we cannot yet defend.

## 2. Close the BP8 convergence question

Edge/centre slip ratio moves with domain extent far more than with resolution,
and straddles the Kim (HBI) reference rather than converging on it:

| x,z extent | y extent | slip(0,0) | slip(−200) | edge/ctr |
|---|---|---|---|---|
| ±500  | ±500  | 37.31 | 21.99 | 0.589 |
| ±1500 | ±500  | 36.90 | 19.52 | 0.529 |
| ±500  | ±2000 | 38.33 | 29.66 | 0.774 |
| Kim   | —     | ~38   | ~21   | 0.553 |

Two sweeps in flight: x/z extent at ±1000 and ±3000, and the time-step factor
`xi` at 0.4 / 0.1 / 0.05 against the current 0.2. Note `xi` sets
Δt = ξ·D_RS/V (`faulting.f90`), not a CFL condition — cell size does not enter,
so refining `dx` does not shrink the step.

Watch for: `fxmin/fxmax` and `fzmin/fzmax` set **both** the elastic domain and
the fault extent. At ±500 the fault cuts edge-to-edge through the model and the
boundary sits 100 m outside Ω_f (800 × 800 m). Widening x/z adds locked fault,
not new frictional area, which is the intended approximation to an unbounded
medium — but it is an approximation, and it is the thing being tested.

If ±3000 lands near ±1500, we have a converged answer and the 10 m run is
worth it. If it keeps drifting, refining is wasted and something structural
needs finding first.

## 3. Multi-fault in the gate set — DONE (v1.11.0/v1.12.0)

`bp1002.qdc.2500` has a reference (`reference/bp1002/cycle0/`), a row in
`tests/e2e/cases.py`, and a physical invariant in
`tests/unit/test_physical_invariants.py` asserting the seed lands on fault 0
and only fault 0. `PROJECT_RULES.md` rule 12 cites a check that exists.

It sits in the **full** tier, not every-push CI: 3821 steps, ~2600 s on 3
ranks against a GitHub runner's four cores. So the fast tier still has no
`ntotft > 1` case. That gap is costed and recorded, not forgotten.

It earned its place immediately. Two defects surfaced only because this case
existed: the fault plane clamped across the whole model
(`meshgen.f90`/`mesh4num.f90` tagging on y alone, 1222 of 1528 nodes per
plane frozen), and two of three requested stations silently matching no node.
Both were invisible to every single-fault benchmark.

## 5. Step-over science — answered, and the answer is conditional

Run and re-run. The result moved three times, which is itself the finding:
the first two answers were artefacts.

1. Rupture appeared **not** to cross. Artefact: the clamped fault plane
   (item 3) manufactured +25 MPa of normal-stress change.
2. With the mesh fixed, rupture **crossed** in cycle 0, both segments slipping
   ~9.4 m. But that geometry left ~51 % of each segment velocity-weakening
   with thin margins.
3. Rebuilt on BP5's box and BP5's zoning (62.5 km segments, VW 1062 km²
   against VS 1412 km², near BP5's own ratio), cycle 0 does **not** cross.

The multicycle sequence is what settles it. Five cycles of the v1.12.0
geometry:

| cycle | interval | peak V | slip A | slip B | |
|---|---|---|---|---|---|
| 0 | — | 0.846 | 4.99 m | 0.03 m | A alone |
| 1 | 0.06 yr | 0.899 | 0.92 m | 4.88 m | B alone |
| 2 | 478 yr | 0.347 | 15.08 m | 15.08 m | **through-going** |
| 3 | 0.003 yr | 0.023 | 1.59 m | 0.00 m | small, A only |
| 4 | 54 yr | 0.371 | 2.23 m | 5.02 m | B-dominated |

**The step-over is a conditional barrier.** It stops a rupture arriving with
cycle 0's stress state and does not stop one arriving after 478 years of
loading. Nothing drawn from cycle 0 alone describes the system.

Still open, and the reason the margin claim is supported but not isolated:
the VW extent is written in global `|x|`, so zoning and segment geometry are
not independent variables — moving the segments changes what the rule does to
them. A clean experiment defines VW relative to each segment.

Also open on this case:

- The interior tips facing across the step-over are velocity-weakening and
  untapered. That drove effective normal stress to zero and STOP 508 by cycle
  3 in the pre-v1.11.0 geometry. It has not recurred in five cycles of the
  current one, but it is an idealisation, not a result.

## 5b. Mid-cycle resume is not possible

`fault.r.nc` carries 12 variables -- shear stress (strike, dip), effective
normal stress, slip rate, state, state-normal, and the six master/slave
velocity components -- but not accumulated slip, and simulated `time` is never
restored: it initialises to 0 and accumulates via `time = time + dtev1`
(`solveTimeLoopMUMPS.f90:72`).

That is correct for its intended use. A *cycle* restart legitimately resets the
clock and the slip accumulator, since output is per-cycle in `Q0/`, `Q1/`,
`Q2/`, and it is the stress and state carrying over that make the next cycle
physical.

But a run that hits `nstep` mid-cycle cannot be resumed. Stress and state
continue correctly while time and slip both reset to zero, so an event ends up
split across two files that cannot be concatenated -- the second file's `t = 0`
is not the first file's end, and its slip starts from nothing. Until this is
fixed, `nstep` must be large enough to reach the slip-rate exit in one run.

Adding `time` and `slips`/`slipd` to the restart file, restored when
`icstart > 1`, would fix it and would also let a long run survive a box reboot.

## 6. Known limits, named rather than hidden

- `porepressure.f90` indexes `nftnd(1)` — no fluid injection on a multi-fault
  model, so BP8-style problems remain single-fault.
- `case.setup` writes one station list to every fault, and
  `library_output.f90:11` writes stations only under `if (j==1)` — faults 2
  and beyond get no station output at all. Not aliasing: they do not exist.
- The post-processing utilities report the **global** peak slip rate, so
  which segment ruptured is invisible in `peak_slip_rate_vs_time` and
  `accumulatedSlip`. On a step-over that is the question being asked.
- `plotOnFaultVars`' depth axis is the fault's own extent with no indication
  of the domain's, which reads as "the model stops at 20 km" for bp1002.
- Reference version skew: bp1002 1.12.0, bp5 1.7.2, bp7 1.7.0, bp8 1.6.0. No
  two are same-binary comparable (rule 11).
- A kink (bent-fault) compset for Liu et al. (2020) is paused in the worktree
  `/home/utig5/dliu/eqquasi.kink`, branch `kink-geometry`, uncommitted. The
  geometry generator is verified in the file and in the solved mesh; dx = 300 m
  (the paper's) needs a 64-bit-integer MUMPS, since the factors want 5.51e9
  reals against a 2^31 ceiling. dx = 600 m runs.
- `func_lib.f90`'s `insert_rough_fault` is single-fault by construction.
- Kim's state variable sits flat at log₁₀θ ≈ 2.8 where ours grows. The aging
  law cannot produce flat, and with the slip law dropped from BP8 on
  2026-08-06 this is unexplained. It is the one question still worth putting to
  Kim — the two spec gaps that were in the draft email (over-determined initial
  condition, missing Poisson's ratio) are both closed by the 2026-08-12
  revision.
- Shared subroutines with `allocatable` or `pointer` dummies need an explicit
  `interface` block in every caller. This codebase uses external subroutines
  with no module wrapper, so an implicit interface silently corrupts the
  descriptor under gfortran and surfaces as a segfault inside the callee.
