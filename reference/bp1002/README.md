# BP1002 — two-fault step-over

BP5 parameters and BP5's frictional zoning, geometry changed only: two
right-lateral segments offset across strike. This is the gate's **only
`ntotft > 1` case**.

| | |
|---|---|
| Segment A | `x ∈ [-60, 2.5]` km at `y = 0` km, cuts the mesh at `x = -60` km |
| Segment B | `x ∈ [-2.5, 60]` km at `y = -5` km, cuts the mesh at `x = +60` km |
| Both | `z ∈ [-20, 0]` km, `dx = 2500` m, 234 nodes each |
| Overlap | 5 km along strike, `x ∈ [-2.5, 2.5]` km |
| Domain | `x ±60`, `y ±50`, `z -60…0` km — BP5's box |

Each segment is 62.5 km long and they overlap by 5 km, so the along-strike
extent of the model is BP5's 120 km and both segments run to a lateral
boundary.

## Zoning, and why it decides the answer

BP5's, expressed in the model's own `x`:

| zone | extent | area |
|---|---|---|
| non-rate-and-state, fixed creep | `\|x\| > 50` km | 450 km² |
| velocity strengthening | `\|x\| ≥ 40` km | 1412 km² |
| velocity weakening | `\|x\| ≤ 38` km and `4 ≤ \|z\| ≤ 16` km | 1062 km² |

with 2 km transitions, as BP5 uses. For scale, BP5 itself has 868 km² of VW
inside 3416 km² of VS — a 4:1 margin — and its rupture arrests after breaking
22% of its fault nodes.

**In cycle 0, the margin decides whether rupture crosses — not the
step-over.** Holding this zoning fixed and shortening the segments to 40 km —
which raises VW from 36% to 49–56% of each segment while leaving the total VW
area almost unchanged, 1000 vs 1062 km² — makes the first rupture cross:
fault B takes 12.03 m on all 153 nodes, peak 2.394 m/s. With 62.5 km segments
it does not.

Read that as a statement about the first event only. The multicycle sequence
below shows the same 62.5 km geometry producing a through-going rupture in
cycle 2, so the margin sets what the *initial* stress state can drive across,
not whether the step-over is passable at all.

A caveat on that comparison: the zoning is written in global `|x|`, so zoning
and segment geometry are not independent variables — moving the segments
changes what the rule does to them. An experiment that separated the two
cleanly would define the VW extent relative to each segment.

## What cycle 0 does

Fault 0 nucleates from the seeded patch, ruptures its velocity-weakening
region, and arrests inside its own segment — the outer ~20 km never slips.
Fault 1 does not rupture: 0.03 m on 10 of its 234 nodes, creep-level.

| | fault 0 | fault 1 |
|---|---|---|
| max slip | 4.99 m | 0.03 m |
| nodes > 1 cm | 169/234 | 10/234 |
| peak slip rate (global) | 0.846 m/s | |
| Δσ̄ | +2.97 … −0.06 MPa | |

3821 steps, ~2600 s on 3 ranks, 41472 elements.

**Cycle 0 is not representative, and a single-cycle reference cannot show
it.** A 5-cycle run of the same compset gives a mixed sequence:

| cycle | interval | peak V (m/s) | slip A | slip B | |
|---|---|---|---|---|---|
| 0 | — | 0.846 | 4.99 m | 0.03 m | A alone |
| 1 | 0.06 yr | 0.899 | 0.92 m | 4.88 m | B alone |
| 2 | **478 yr** | 0.347 | **15.08 m** | **15.08 m** | **both — through-going** |
| 3 | 0.003 yr | 0.023 | 1.59 m | 0.00 m | small, A only |
| 4 | 54 yr | 0.371 | 2.23 m | 5.02 m | B-dominated |

Cycle 2 crosses: after 478 years of loading both segments slip 15.08 m,
equal to two decimals, on 196 and 194 of their 234 nodes. That is one rupture
spanning both faults, not two events.

So the step-over is a conditional barrier, not an absolute one. It stops a
rupture that arrives with cycle 0's stress state and does not stop one that
arrives after a long interseismic period. Any statement drawn from cycle 0
alone -- including the margin comparison above, which is also a cycle-0
comparison -- describes the first event from an artificial initial condition,
not the system's behaviour.

## Why this case is in the gate

Every other benchmark has `ntotft = 1`, so nothing else exercises per-fault
routing of on-fault input, and three multi-fault bugs in this project's
history were found by hand for that reason (`PROJECT_RULES.md`, rule 12).
`case.setup` writes a `(ntotft, nfz, nfx, nvar)` array,
`netcdf_read_on_fault` reads it back with an explicit fault dimension, and
`testsys/unit/test_physical_invariants.py` asserts the seed lands on fault 0
and only fault 0.

## Two things that bit, recorded so they do not again

**The fault plane was clamped.** `meshgen.f90`/`mesh4num.f90` tagged every
node on a fault's y-plane as a fault boundary on the y coordinate alone, so
1222 of the 1528 nodes on each plane sat frozen at zero displacement while
the surrounding material moved metres. It manufactured +25 MPa of
normal-stress change and inverted the physics. Fixed in v1.11.0; any bp1002
result from before that is superseded.

**Two stations vanished.** An earlier version of this compset asked for
stations at `x = ±2.0` km and `z = -2.0/-18.0` km, none of which is a
multiple of `dx = 2500` m. They matched no node and were dropped with no
message, leaving one station where three were requested — and the reference
was blessed from it. `src/eqquasi.f90` now reports the mismatch, and the
coordinates here are on-grid.

## Why dx = 2500 m and not BP5's 2000 m

The fault-normal offset is 5000 m and the mesher (`src/func_lib.f90`,
`build_yline_belt`) requires it to be an integer multiple of `dy = dx`:
`5000/2000 = 2.5` is refused, `5000/2500 = 2` meshes.

## Tier

`cycle0/` — the complete event, **full** tier only. Not in every-push CI:
3821 steps and ~2600 s on 3 ranks of a 64-core host, against a GitHub
runner's four cores. That leaves the fast tier with no `ntotft > 1` case; the
gap is costed and recorded in `testsys/e2e/cases.py`, not an oversight.

**Caveat.** First reference for this benchmark — a regression lock, not a
validation. It has not been checked against any independent oracle, and the
interior tips facing across the step-over are left velocity-weakening and
untapered, an idealization that drove the effective normal stress to zero by
cycle 3 in an earlier geometry.

## Provenance

Blessed at **v1.11.0** (`cycle0/runInfo.json`). Confirmed still reproducible
by v1.13.0 on 2026-08-14: all 19 pre-existing entries compared equal under
the e2e comparators, no exceptions.

The `fltst_ft2_*` files were added on 2026-08-14 and are the only part of
this reference not produced at v1.11.0. Until then `output_onfault_st` opened
a station file only for fault 1, so the three stations this case requests on
fault 2 wrote to an unnamed `fort.51` and never reached the reference at all;
the gate walks reference files, so their absence made them silently
uncompared rather than visibly missing. They were added additively, and only
after the fault-1 files were confirmed byte-for-byte unchanged, so nothing
here was re-blessed to accommodate a code change.

In cycle 0 fault 2 barely moves -- these stations record 1e-6 to 1e-4 m of
slip creeping at ~1e-9 m/s, against fault 1's 1.65 m. That is the point of
keeping them: the step-over holding is a result, and until now the gate had
no record of the segment that did not rupture.

## How entangled the zoning and the geometry are

Zoning is written in the model's global |x|, so a segment's frictional
makeup depends on where it sits, not on what it is. Measured on the nodes
each run actually meshed:

| run / segment | VW | VS | non-RSF |
|---|---|---|---|
| current, A (−60 → +2.5 km) | 36.3% | 48.3% | 15.4% |
| current, B (−2.5 → +60 km) | 36.3% | 48.3% | 15.4% |
| `ab_newzone_olddom`, A (−42.5 → −2.5 km) | 49.0% | 51.0% | 0% |
| `ab_newzone_olddom`, B (−7.5 → +32.5 km) | 55.6% | 44.4% | 0% |

Two things follow, and the second was not previously stated.

The current geometry is symmetric about x = 0, so both segments receive
*identical* zoning. Whatever decides which one ruptures, it is not a zoning
difference between them.

The A/B experiment is a different matter. It was read as isolating segment
length, and it does not: shortening the segments also raised VW from 36% to
49–56% and removed the non-rate-and-state creeping ends entirely. Worse, its
two segments are not zoned alike — B carries 6.6 points more VW than A —
because that geometry is not symmetric about the step-over. So "the VS margin
decides cycle-0 crossing" remains supported but not isolated, and the
asymmetry means which fault went first in that run was partly a zoning
choice.

Making it isolated needs VW defined in each segment's own along-strike
coordinate, holding VW fraction and non-RSF fraction fixed while segment
length varies. That is a new compset, not a parameter change.

## The 9-cycle sequence, and how it ends (2026-08-15, v1.14.0 binary)

A 20-cycle run at `nstep = 30000` (`work/bp1002_stepover/multicycle_20_knox`)
ended at cycle 9 on the non-compressive-normal-stress guard (`STOP 508`):
sigma_bar reached +0.60 MPa at x = −2.5 km, z = 0 — fault A's interior tip at
the free surface, facing the step-over — 87 yr into the cycle. Every one of
the nine completed cycles ended on its own physics, none on the step cap.

| cycle | steps | peak V (m/s) | slip A / B (m) |
|---|---|---|---|
| 0 | 3821 | 0.846 | 4.99 / 0.03 |
| 1 | 4572 | 0.899 | 0.92 / 4.88 |
| 2 | 11171 | 0.347 | **15.08 / 15.08** |
| 3 | 5854 | 0.375 | 2.23 / 5.03 |
| 4 | 11362 | 0.739 | **12.05 / 12.05** |
| 5 | 6627 | 0.817 | 2.73 / 5.94 |
| 6 | 10885 | 0.338 | **12.56 / 12.56** |
| 7 | 8567 | 0.776 | 6.24 / 6.90 |
| 8 | 8542 | 0.785 | **7.34 / 7.34** |

Two observations the truncated 10000-step baseline could not show:

- **Through-going events alternate with confined ones** (bold rows), and their
  slip declines monotonically: 15.08 → 12.05 → 12.56 → 7.34 m.
- **The run's end is the caveat above coming due.** The untapered
  velocity-weakening interior tips accumulate a normal-stress reduction with
  every through-going rupture; four of them progressively unclamped the tip
  until rate-and-state had no solution. The guard stopped the run rather than
  silently producing NaNs. Whether to taper the tips — which changes the
  model — is an open scientific decision, not a bug fix.

The old `nstep = 10000` baseline truncated cycle 2 mid-descent, so its
cycles 3+ (restarted mid-event) are artifacts; this sequence supersedes it
from cycle 2 onward. Cycles 0–1 match it to 1e-12.
