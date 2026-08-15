# Pathway forward

Working plan, ordered by what blocks what. Updated 2026-08-14, at v1.13.0.

This is a live document: close items by deleting them, and record what was
learned in the place that will be read again (`PROJECT_RULES.md` for a rule,
`reference/<bench>/README.md` for a benchmark finding, a commit message for a
fix). Do not let this file become the archive.

---

## 55-HOUR AUTONOMOUS QUEUE — started 2026-08-14 17:00, user back Monday

Read this first on wake-up. Update in place; close items by deleting them.

### Standing rules
- **Read `PROJECT_RULES.md`'s action card before staging and before adding any
  file, name or parameter.** Every rule broken here was already written when it
  was broken. The card is 20 lines at the top of that file.
- Host is **knox** (64 cores). Check `uptime`; at most 2-3 runs at 3 ranks;
  wait above load 56. cotopaxi is a different machine and unreachable from here.
- **Never rebuild while a multi-cycle run is in flight** (rule 9). run.sh
  invokes `eqquasi-<version>` afresh each cycle, so a mid-run rebuild splits
  the dataset. Wait for a gap, or bump the version so the filename differs.
- Background commands inherit a stale cwd. Always `cd` to an absolute path
  inside the command.
- Never `pkill -f <pattern>` matching the issuing shell. Kill by PID via
  `readlink /proc/<pid>/cwd`.
- `gh` is not on PATH: use `/home/staff/dliu/bin/gh`.
- Two runs of the same binary are never byte-identical (MPI reduction order,
  1e-15 to 1e-21). Compare numerically, never by hash.

### Queue

1. [x] **Release: v1.15.1, CI GREEN (run 31858202391).** Recorded honestly:
   victor cut v1.15.0 on a master that was already red (rule 17 breach -- the
   workflow PATH still exported the deleted scripts/), and his audit found the
   caps read violated rule 4: on a pre-caps model.txt with faultgeom it
   consumed fault 1's geometry line (reference/bp1002's file parsed to
   min_norm=-60000 Pa). v1.15.1 fixes the PATH, reads the caps record whole
   (backspace when it is faultgeom, loud stop 8 when malformed), runs
   checkNormalStressCaps every cycle instead of cold start only, and stops on
   missing restart inputs. Verified live both ways. All four live runs carry
   the caps line; eqquasi-1.14.0 untouched. v1.15.0's tag stands, superseded.

2. [~] **Runs.** bp1002_stepover/multicycle_20_knox ENDED AT 9 OF 20 CYCLES
   on the stop-508 physics guard: effective normal stress went non-compressive
   (+0.60 MPa) at x = -2.5 km, z = 0 -- fault A's interior tip at the free
   surface, facing the step-over -- 87 yr into cycle 9. This is the
   documented idealization (untapered VW inner tips); four through-going
   events (c2 15.08 m, c4 12.05, c6 12.56, c8 7.34, declining) progressively
   unclamped the tip. The 9-cycle dataset STANDS; fixing it means tapering
   the tips, a model change for the user. Kink dtcap verdict: 5 cycles,
   1369 yr, flat 1.173e-9 through the epoch where uncapped dtmax=0 spiked to
   98 m/s -- uncapped adaptive dt manufactured those spikes; a default cap is
   a solver-behaviour change, USER'S CALL. dip90 control: bend exonerated
   (3 physical cycles). kink600: c0-c1 coherent, 0 swings, but every cycle
   hits the nstep=10000 cap; gracefully stopped after the in-flight cycle to
   free the box for the e2e tier.
   - `work/bp1002_stepover/multicycle_20_knox` — 20 cycles, nstep=30000.
     Cycles 0-1 match the 1.13.0 baseline exactly; cycle 2 ran 11171 steps and
     exited on physics at 7.2e-6, where the old nstep=10000 run truncated it.
     This is the science run. Compare each cycle with
     `work/bp1002_stepover/compare_cycles.py multicycle_20 <run>`.
   - `work/liu2020.kink.10cyc` — dx=600, Lambda0/dx=0.94. Produced a coherent
     rupture (rise to 2.2e-1, then decay) where dx=1200 oscillated. Watch
     whether later cycles stay coherent or degrade as 1200 did.
   - `work/kink.dc014.dx1200` — Dc=0.14, Lambda0/dx=6.0. Controlled test: same
     mesh and geometry as the oscillating run, one parameter changed. Still
     interseismic at 1.7e-9 after 102 blocks, as expected for the larger
     nucleation size.
   - `work/liu2020.kink.10cyc.dx1200` — finished, 10 cycles, oscillating.
     Kept as the negative control. Do not delete.

3. [ ] **Full e2e tier at HEAD.** Has not run since the compsets rename, the
   script/testsys rename, or the cohesive-zone precheck. `MACHINE=utig` on
   cotopaxi, `MACHINE=ubuntu` on knox -- knox has no `utig` branch in
   install.eqquasi.sh and falls through to no LD_LIBRARY_PATH. ~3h05.

4. [x] **Five pre-par compsets ported** (values untouched, dy/dz paired --
   the contract test caught the old files never set them). UNVERIFIED: no
   oracle. Surfaced+fixed: generateFaultInterface wrote geometry to the case
   root while the solver reads input/ (latent since v1.13, ungated because
   test.bp5.qdc.dip90 is orphaned); rough tables renamed to the read name.
   Was: **The five broken compsets.** bp1001.fdc.250, bp1001.fdc.rough.250,
   bp1001.qdc.rough.250, liu2020.fdc.planar.300, liu2020.fdc.rough.250 declare
   bare module-level variables instead of the `par` object, so `case.setup`
   dies on `NameError: name 'par' is not defined`. Porting them is safe to
   attempt ONLY as far as making them run; none has a reference, so a port
   cannot be verified and must not be blessed. If a port runs, say explicitly
   that it is unverified.

5. [ ] **fric() index revamp** — named constants in one place, applied
   everywhere, plus an index -> meaning -> unit -> writer table. Must be one
   pass: naming at some call sites only creates a fourth source of truth on
   top of rule 5's three. Rule 5 needs rewriting when it lands.

6. [ ] **Per-fault peak slip rate.** `global.dat` column 2 is a single maxval
   over all faults in solveTimeLoopMUMPS.f90. Per-fault peaks need new columns
   and a re-bless of every reference. DO NOT DO THIS UNATTENDED -- it is an
   output-format change the user must approve.

7. [ ] **Orthogonal zoning experiment** — needs a new compset (rule 1: flag
   first), and holding VW fraction fixed vs VW area fixed are two different
   experiments answering different questions. User's call, not a default.

### Known open, not queued
- GitHub Releases stop at v1.3.2 (2022) while 32 tags exist. Tags are not
  Releases; cutting them going forward is easy, backfilling 26 means inventing
  notes.
- `src/globalvar.f90` has one comment reading `scripts/case.setup`, left stale
  because correcting it would invalidate the binary three live runs are using.
- 243 `.DS_Store` files across `$HOME`, none in this repo any more.
- `test.bp5.qdc.dip90` and `test.stepover.qdc` are regression variants with no
  e2e row -- orphaned.

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
- **Hardcoded normal-stress caps in `src/faulting.f90:169-177`.** Known, not
  being fixed now.

      max_norm = -40.0d6
      min_norm = -10.0d6

  active whenever `rough_fault == 1 .and. C_elastic == 1`. Liu, Duan & Luo
  (2020) §3.5 states the caps as **-100 to -10 MPa**; the -10 matches, the
  -40 does not. At -40 the solver clamps that paper's own -50 MPa initial
  normal stress on the first step, and Fig. 6's excursion to -100 MPa at the
  bend cannot be represented at all.

  Neither number is a parameter -- both are literals, in no compset and no
  docstring. Five compsets reach this branch (bp1001.fdc.250,
  bp1001.fdc.rough.250, bp1001.qdc.rough.250, liu2020.fdc.rough.250,
  liu2020.kink.qdc), so the values cannot simply be changed: they need
  parameterising with -10/-40 as the default so existing cases are untouched,
  and -10/-100 set for the kink compset. Hardcoded physics thresholds that
  silently override a compset's stated initial condition are the wider
  problem; this is one instance.

- **Reproducing Liu, Duan & Luo (2020), EQquasi half.** Paper at
  `work/liu2020.kink/paper/Liu_Duan_Luo_2020_GJI.pdf`. Read pp. 1-8; the
  appendix and switch criterion (pp. 9-12) matter only for the coupled EQsimu
  loop, not the EQquasi-only compset.

  Table 1: FS 60 km, FD 30 km, VW ~45 km long x 11.4 km wide, V_min 1e-12,
  Vs 3464, Vp 6000, rho 2670, V0 1e-6, mu0 0.6, b 0.011, a 0.007, sigma_n
  -50 MPa, tau 30 MPa, L 11 mm, dx 300 m, dt_min 0.025 s.

  §3.5 initial stress: uniform regional field sigma_1 = -80 MPa, sigma_3 =
  -20 MPa at 45 deg to the left segment, giving shear 30 MPa / normal -50 MPa
  there; the right segment, rotated 10 deg, gets normal -60.26 MPa and shear
  set to 33 MPa for consistency. Starts from uniform slip rate 1e-9 m/s at
  steady state.

  §3.3 friction: a, b depth-dependent multilinear per Lapusta & Liu (2009);
  b reduced where |x| > 23 km and where |z| > 17.5 km or |z| < 4 km; L = 11 mm
  rising linearly to 27 mm between 4 km depth and the free surface.

  §3.2 EQquasi domain: x, y -30..30 km, z -30..0 km (NOT EQdyna's 120x120x60);
  y boundaries driven at +/-9.513e-10 m/s giving ~1e-9 m/s creep; other faces
  free; dx 300 m uniform in x, z, enlarged 1.3x in y beyond 1.2 km.

  **§3.4 resolution.** The paper requires `h*/W_seis < 1` and
  `Lambda_0/dx > 2.3`, and reports 0.8674 and 2.5163 at dx = 300 m.
  Lambda_0/dx scales as 1/dx, so 600 m gives 1.26 and 1200 m gives 0.63 --
  both below the criterion. Real, and it would matter for a rupture front.

  **It is NOT what made the smoke test blow up.** An earlier version of this
  note blamed resolution; that was wrong. The trace says otherwise: step 1 is
  correct at 1.17e-9 m/s, step 2 jumps 14 orders of magnitude in moment rate,
  and time then freezes -- a one-step instability from a steady creeping
  initial condition, which under-resolution does not produce.

  The caps do. `user_defined_params.py` computes initial shear from
  `shear_steady_state()` at sigma_n = -50 MPa; `faulting.f90` then clamps
  sigma_n to -40 MPa on the first evaluation and never recomputes shear. The
  first snapshot shows effective normal stress spanning exactly 40.00 to
  50.00 MPa and tau/sigma reaching 0.785 against mu0 = 0.6. With a = 0.007,
  RSF turns that into V/V_ss = exp(0.15/0.007) ~ 2e9, which takes 1e-9 m/s to
  ~2 m/s in one step -- exactly what step 2 did. The caps do not merely
  disagree with the paper; they destroy the compset's own initial condition.

  BLOCKER: dx = 300 m is not optional, and there the MUMPS factors want
  ~5.51e9 reals against a 2^31 ceiling, so full reproduction needs a
  64-bit-integer MUMPS build.
- **fric() index revamp — systematic naming and documentation.** `fric(1:100,
  node, fault)` is indexed by bare integers everywhere: `fric(26,...)` is
  trial slip rate, `fric(46,...)` V_init, `fric(81:84,...)` the rupture-area
  block, and nothing in the source says so except trailing comments that only
  some sites carry. Rule 5 currently *documents* this as the state of affairs
  and guards additions (pick an unused index, write it in three places); it
  does not fix it.

  The revamp: named constants in one place, applied everywhere, plus a table
  mapping index -> meaning -> unit -> which of the three files writes it.
  It must be done in one pass. Introducing names at some call sites only
  creates a fourth source of truth on top of the three rule 5 already names
  (`script/defaultParameters.py`, `script/case.setup`, `src/netcdf_io.f90`),
  which is worse than the bare integers. Rule 5 will need rewriting when it
  lands, since it exists to describe the unnamed state.
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
4.3 (`script/resampleBP8Profiles.py`), validate (`script/checkBP8Submission`),
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
`testsys/e2e/cases.py`, and a physical invariant in
`testsys/unit/test_physical_invariants.py` asserting the seed lands on fault 0
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
