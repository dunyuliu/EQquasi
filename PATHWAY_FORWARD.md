# Pathway forward

Working plan, ordered by what blocks what. Updated 2026-08-14, at v1.13.0.

This is a live document: close items by deleting them, and record what was
learned in the place that will be read again (`PROJECT_RULES.md` for a rule,
`reference/<bench>/README.md` for a benchmark finding, a commit message for a
fix). Do not let this file become the archive.

---

## 24-HOUR AUTONOMOUS QUEUE — started 2026-08-14, unattended

Read this first on wake-up. Update in place; close items by deleting them.

### Standing rules
- Verify the binary before every run: `grep EQQUASI_VERSION src/globalvar.f90`
  against `strings bin/eqquasi-<v> | grep Welcome`. Rebuild with
  `MACHINE=utig`, install as `bin/eqquasi-<version>` (no plain `eqquasi`).
- Background commands inherit a stale cwd. ALWAYS `cd` to an absolute path
  inside the command; this has broken four launches today.
- Never `pkill -f <pattern>` where the pattern matches the running command.
- Killed runs leave NFS `.nfs*` files; kill the ranks, sleep, then `rm -rf`.
- Do not work in the user's home outside this repo. The kink work lives in
  the worktree `/home/utig5/dliu/eqquasi.kink`, NOT in `~`.
- No new files or names without need (rule 1). No force-push.
- **Respect the shared host (rule 15).** This is a 64-core machine other
  people are using: baseline load was ~50 at handoff with only 6 of those
  ranks mine. So:
    * at most TWO of my runs at once, 3 ranks each -- never more;
    * check `uptime` before launching anything, and if load is above ~56,
      wait rather than add to it;
    * never raise HPC_ncpu above 3 for a queue run;
    * a long run is not more urgent than someone else's interactive session.
  Serialise the heavy items (5-cycle, e2e full tier) rather than running them
  together.

### Queue
1. [x] **5-cycle bp1002 on the new workflow reproduces the science.**
   All five cycles run, same event sequence, same slips:

   | cycle | steps | peak V | rel diff | slip A | slip B | slip rel diff |
   |---|---|---|---|---|---|---|
   | 0 | 3821/3821 | 0.845757 | 2e-15 | 4.99 m | 0.028 m | 0 |
   | 1 | 4572/4572 | 0.899183 | 8.7e-12 | 0.919 m | 4.876 m | 0 |
   | 2 | 10000/10000 | 0.347357 | 7.7e-9 | 15.08 m | 15.08 m | 0 |
   | 3 | 984/984 | 0.022665 | 2.7e-7 | 1.587 m | 0.0002 m | 6.5e-8 |
   | 4 | 5845/**5846** | 0.370585 | 2.4e-6 | 2.226 m | 5.023 m | 9.5e-8 |

   Agreement degrades by roughly a factor of 1000 per cycle, and by cycle 4
   the run exits one step earlier. That is accumulation, not a defect, and
   cycle 0 is what proves it: identical inputs, identical binary, only the
   file paths changed, and it still differs by 2e-15. A single cycle is
   already nondeterministic at round-off because MPI reduction ordering is
   not fixed, and each cycle starts from the previous cycle's restart file,
   so the seed is amplified by a system with a positive Lyapunov exponent.

   The signature discriminates the two explanations. Accumulated round-off
   grows smoothly from 1e-15; a broken restart handoff -- a dropped variable,
   the wrong file, a stale path -- would show as a step change at cycle 1,
   the first restart, and stay there. Cycle 1 agrees to 8.7e-12. Nothing
   jumps anywhere.

   Do not read this as bit-reproducibility, which this code does not have and
   would need a fixed reduction order to get. Read it as: the workflow revamp
   moved no physics.

   Compare with `work/bp1002_stepover/compare_cycles.py <baseline> <run>`.
   Two traps it now avoids, both of which cost time here:
   - Peak V is column 2 of global.dat (every step). The max over
     fault.NNNNN.nc snapshots understates it ~20x -- they are written every
     1000 steps and step over the peak.
   - The baseline is `multicycle_20`, NOT `old_5cyc_prefix`, which predates
     the meshgen fault-plane fix and gives completely different numbers.
     Comparing against it reads as a total regression when nothing is wrong.

2. [x] **e2e gate green, both tiers.** Fast: 23 passed, 6 skipped, 0 failed.
   Full: 41 passed, 13 skipped, 2 failed -- both diagnosed and fixed, and the
   fast tier re-run green afterwards. The full tier takes 3h05 and needs
   `MACHINE=utig` on this host.

   The two failures were both gate defects, not physics:

   a. `bp5 offfault srfst_strk000st032dp000.txt`: 1 entry of 31381 differing
      in the 7th printed digit. The station files are written `E21.13,6E15.7`
      -- seven significant figures -- while the gate demanded rtol=1e-9. That
      asks the file for digits it never wrote, so MPI reduction-order noise
      (~2e-14) flipping one rounding boundary failed the suite at random
      instead of on regressions. `column_printed_rtol()` in tests/e2e/cases.py
      now floors each column's rtol at 2 ulp of its own printed precision,
      read off the file rather than hand-tuned. Verified: the real run passes
      at 9.9e-17 and a 1e-5 relative perturbation is still caught. Anything
      below the 7th digit is undetectable in this tier by construction -- the
      netCDF files carry full double precision and are still compared at 1e-9.

   b. `test_clean_build.py` pointed at `bin/eqquasi`, which stopped existing
      when binaries became versioned. So `test_built_binary_version_matches_
      the_source` -- the stale-binary guard added *because* of the 1.7.0 vs
      1.10.0 incident -- skipped itself with "no binary" on every run while
      the suite reported green. It now derives `bin/eqquasi-<declared>` and
      treats a missing binary as a failure, not a skip.

   Earlier note here blamed (b) on "no system MUMPS on this host". Wrong: the
   test defaults MACHINE to ubuntu, and the ubuntu target looks for
   dmumps_struc.h where this host does not keep it. The build works fine with
   MACHINE=utig. The assertion now names MACHINE so the next reader is not
   sent looking for a missing library.

3. [x] **Revamp committed, released as v1.13.0, pushed, CI green.**
   Six commits, tagged v1.13.0. CI run 31785517943 SUCCEEDED -- the first
   green since v1.12.0.

   CI had been red across every push since v1.12.0, and the note that used to
   sit here calling the only failure "local-only, no MUMPS" was wrong twice:
   wrong about the local cause, and wrong that remote CI was fine. The
   workflow's Build step asserted `test -x bin/eqquasi`, the unversioned name
   that stopped existing when binaries became versioned, so it died at Build
   before the e2e tier ran at all.

   That is three places one change broke, all found the same day: this, the
   e2e stale-binary guard, and the contract test's filename split. When a
   name that appears in build products changes, grep the whole tree for the
   old one -- tests and CI config included, not just source.

4. [~] **Then work the list below**, cheapest first. Do not start anything in
   section 1 or 2 (BP8 submission/convergence) unattended -- they need
   judgement calls on what to submit.

5. [x] **Kink work staged, not merged.** Committed as 7d11983 on
   `kink-geometry` in the worktree `/home/utig5/dliu/eqquasi.kink`. Not
   pushed, not merged.

   Read before picking it up: that branch sits at e8327a7 and PREDATES the
   v1.13.0 workflow revamp, so the compset has never been run against a case
   with input/, result/, scratch/. It was staged rather than brought forward
   because rebasing it onto master and re-verifying is real work with a real
   chance of breaking, and the instruction was not to merge unattended.
   Geometry is verified in the written file and in the solved mesh; dx=600 m
   runs, dx=300 m (the paper's) needs a 64-bit-integer MUMPS -- the factors
   want 5.51e9 reals against a 2^31 ceiling.

### Landed during the revamp, worth remembering
- **Implicit interfaces bite character arguments too.** Passing
  `trim(inputDir)//"name.nc"` to an external subroutine whose dummy is
  `character(len=50)` hands over a temporary of the expression's own length;
  with no interface block the callee still reads 50 bytes and appends 22 bytes
  of adjacent memory to the path. Symptom was a netCDF "No such file or
  directory" on `input/on_fault_vars_input.nca  0 <junk>`. Assign into a
  fixed-length local first. The same file already warned about this for
  allocatable dummies -- it applies to character dummies as well.

### Known open, from today
- [FIXED] Multi-fault station output. It was worse than "faults 2+ get
  none": the unopened unit 51 was still *written to*, so Fortran connected it
  to a default `fort.51` and every fault-2+ station appended there,
  interleaved. BP1002 asks for six stations and three produced nothing
  usable. Fault 1 keeps the plain SEAS name; faults 2+ take an `ft<N>_` tag.
  CONSEQUENCE CLOSED 2026-08-14: `reference/bp1002` now carries the three
  `fltst_ft2_*` files. Added additively after all 19 pre-existing entries
  compared equal at v1.13.0, so nothing was re-blessed to accommodate the
  fix.
- Utilities report the GLOBAL peak in `peak_slip_rate_vs_time`, so which
  segment ruptured is invisible. NOT a plotting fix: `global.dat` column 2 is
  a single `maxval` over all faults, computed in `solveTimeLoopMUMPS.f90`.
  Per-fault peaks mean writing per-fault columns, and columns 5-7 of
  global.dat are currently zeros that every benchmark reference has recorded
  -- filling them fails the e2e comparison for bp5, bp7 and bp1002 until all
  of them are re-blessed. That is a judgement call on an output format, so it
  is left for the user rather than done unattended.
  (The `accumulatedSlip` half of this item was stale: `plotAccumulated`
  already takes `--fault`.)
- [FIXED] `plotOnFaultVars`' depth axis now names the domain floor when the
  fault does not reach it.
- [FIXED] Reference provenance. Each README now records the version that
  produced its numbers (bp1002 1.11.0, bp5 1.7.2, bp7 1.7.0, bp8 1.6.0) and
  the date current code was confirmed to reproduce them. The version gap was
  never staleness -- the full tier passes against all of them -- it was the
  absence of a record, which is what makes a future divergence datable.
- Every pre-existing case in `work/` can no longer RESUME: the solver
  hard-stops without `input/` and `scratch/`. Deliberate, and narrower than
  it sounds -- verified 2026-08-14 that the utilities still READ the old flat
  layout, `plotPeakSliprateTime.py cycle0 cycle1` and `plotAccumulated
  --fault 1` both included, because cycle discovery tries `result/cycle*`,
  then `cycle*`, then `Q*`. So old results stay analysable; only continuing a
  run needs the new layout.

  No migration script, on purpose: recreating a case is `create.newcase` +
  `case.setup` + copying the last cycle's `*.r.nc` pair into `scratch/` with
  the right `currentcycle.txt`, and a script for that would be a new file
  (rule 1) standing between the user and three commands they can check.
- The orthogonal zoning experiment (VW per segment, not global |x|) is still
  the thing that would make the step-over margin claim isolated rather than
  merely supported. QUANTIFIED 2026-08-14 -- see reference/bp1002/README.md.
  The current geometry is symmetric about x = 0, so both segments get
  identical zoning (36.3% VW each); whatever selects the rupturing segment is
  not a zoning difference between them. But the A/B experiment confounds
  three variables at once: shortening 62.5 -> 40 km also took VW from 36% to
  49-56% and removed the non-RSF ends (15.4% -> 0%), and its two segments are
  not even zoned alike (A 49.0%, B 55.6%) because that geometry is not
  symmetric about the step-over.

  NOT STARTED: it needs a new compset (rule 1 -- flag before creating), and
  the design decision is what to hold fixed. Holding VW *fraction* fixed and
  holding VW *area* fixed are different experiments and answer different
  questions; that is the user's call, not a default.

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
- `case.setup` writes one station list to every fault, so a multi-fault case
  gets the same station coordinates on each. Fault 2+ stations DO exist now
  (v1.13.0); they are named `fltst_ft<N>_*` because the coordinates are
  fault-local and would otherwise collide.
- The post-processing utilities report the **global** peak slip rate, so
  which segment ruptured is invisible in `peak_slip_rate_vs_time`. On a
  step-over that is the question being asked. Not fixable in the plotting:
  `global.dat` column 2 is a single `maxval` over all faults
  (`solveTimeLoopMUMPS.f90`), so it needs new columns and a re-bless of every
  reference. (`plotAccumulated` already takes `--fault`.)
- ~~`plotOnFaultVars`' depth axis gives no indication of the domain's~~
  FIXED v1.13.0: the label names the domain floor when the fault stops short.
- The references were blessed at four different versions (bp1002 1.11.0,
  bp5 1.7.2, bp7 1.7.0, bp8 1.6.0), so no two are same-binary comparable
  (rule 11). Each README now records its version and the date current code
  was confirmed to reproduce it; the gap was never staleness, since the full
  tier passes against all of them.
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
