# Pathway forward

Working plan, ordered by what blocks what. Updated 2026-08-12, at v1.7.0.

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

## 3. Put a multi-fault case in the gate set — highest leverage item here

`BP5`, `BP5.dip90`, `BP7` and `BP8` all have `ntotft = 1`. They prove the
degenerate path is unchanged and say nothing about two faults. Meanwhile
`case_input/bp1002.qdc.2500` and `case_input/test.stepover.qdc` both exist,
both have `ntotft = 2`, and neither is wired into anything.

Order matters, so that `master` never goes red:

1. Fix the violations Zofia's checks detect — wire a multi-fault case into the
   regression/e2e gate, guard `nonfs(i) == 0` in `read_input.f90` /
   `mesh4num.f90` (same silent-degradation shape as the zero-node fault), and
   let the gold agent's unreferenced CSVs resolve as it finishes.
2. *Then* commit `PROJECT_RULES.md` and the three checks
   (`test_reference_gold_is_referenced.py`, `test_gate_set_has_multifault.py`,
   `test_code_convention_landmines.py`).

## 4. The x and z belts carry the bug we just fixed in y

`fltx1 = minval(fltxyz(1,1,:))` is the x belt origin, and other faults' `xlo`
are never checked for commensurability with `dx`. Same in z. v1.7.0 guards only
the fault-normal direction. `bp1002.qdc.2500` happens to be x-commensurate, so
nothing exercises it today. Same refuse-guard, small change.

## 5. Step-over science — the question the case exists to ask

Run `bp1002.qdc.2500` for its three cycles and see whether rupture jumps the
5 km step. BP5 parameters throughout; only the geometry differs, plus a
one-sided velocity-strengthening taper (the interior step-over tips are left
velocity-weakening, an idealisation the owner accepted). The nucleation patch —
BP5's low `Dc` and 0.03 m/s seed — is on fault A's x− end only, so fault B is
unseeded and whether it ruptures is the result, not the setup.

Everything before this has been plumbing.

## 6. Known limits, named rather than hidden

- `porepressure.f90` indexes `nftnd(1)` — no fluid injection on a multi-fault
  model, so BP8-style problems remain single-fault.
- `case.setup` writes one station list to every fault.
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
