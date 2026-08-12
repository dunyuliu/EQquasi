# BP8-QD-GS gold

What this is, what it does and does not establish, and how to regenerate it.

## Configuration

BP8-QD-GS (Gaussian source), aging law, from `case_input/test.bp8.qdc`:

| | |
|---|---|
| cell size | `dx = 50 m` (the benchmark specifies 10 m; §6 invites two resolutions) |
| domain | x, y, z all ±500 m |
| time-step factor | `xi = 0.2` |
| run | 5301 steps to 30.00 days, exits on `fluid_tend` |
| cost | 796 s, 1 rank, 9702 nodes |
| code | v1.7.0 |

This is the configuration packaged as `dliu_eqquasi-1-7-0-50m.zip` and uploaded
to the CRESCENT platform. Gold is that run itself at full time resolution, not a
decimation of it.

## Files

| File | Contents |
|---|---|
| `fltst_strk*.csv` | 9 stations, 5302 rows × 11 columns (§4.1) |
| `global.csv` | 5302 rows × 3 (§4.2: time, log₁₀ max slip rate, moment rate) |
| `*_strike.csv`, `*_depth.csv` | 10 profiles, 664 rows × 19 (§4.3, at the native 17-node grid) |
| `fault.05301.nc` | fault-plane snapshot at the final step |
| `runInfo.json` | node count, steps, wall time, version |
| `summary.json` | per-station scalars the e2e asserts by name, plus provenance |

Written at 17 significant digits. At 9, the time column could not resolve
2.59×10⁶ s to better than 4 ms, which showed up as a uniform 4.167×10⁻³
max\|diff\| across every file — an artifact of the gold writer, not the solver.

## How it is checked

`tests/e2e/test_bp8_against_gold.py` rebuilds, reruns `test.bp8.qdc` at this
configuration, and diffs **every file whole** to 10⁻⁶ — not sampled. Named
scalar assertions sit alongside the full diffs because they say *which* quantity
moved; the full diffs exist because seven scalars out of ~58 000 entries per
station left most of every curve unguarded, including columns nothing sampled
(`slip_3`, both Darcy velocities). The run must also pass
`scripts/checkBP8Submission`.

Quick check without the suite:

```
python3 scripts/plotAgainstGold.py bp8 <run_dir>
```

## What this does not establish

**Gold is a regression lock, not a validation.** It detects unintended change.
It does not establish that BP8 is reproduced correctly, for two reasons:

1. `dx = 50 m` is five times the benchmark's cell size.
2. **The domain is not converged.** x, y, z = ±500 m is the under-converged
   corner. The y faces are rigid walls (`u = 0`) and the fault plane cuts
   edge-to-edge through the box, leaving 100 m outside Ω_f. Edge/centre slip
   moves from 0.589 here to 0.529 once x and z reach ±1000, and to 0.774 when
   the y wall moves to ±2000. See `../README.md` for the full sweep and the
   boundary conditions behind it.

Against Kim's HBI (boundary element, unbounded): centre 37.31 mm vs ~38, edge
21.99 vs ~21, late pressure 1.690 MPa vs ~1.70. Peak pressure 13.63 MPa against
eq. (21)'s analytic 13.0. The agreement at this configuration is partly luck —
the converged numbers are centre 36.9 and edge/centre 0.53.

## Convergence study (independently audited, 2026-08-12)

Sweeps at dx = 50 m, aging law, 30-day runs, under `work/`. Every number
re-derived directly from `fltst_strk*.dat` and `runInfo.json`.

### Time-step factor `xi` — converged

| xi | steps | edge/centre |
|---|---|---|
| 0.4  | 5184 | 0.5893 |
| 0.2  | 5301 | 0.5895 |
| 0.1  | 5927 | 0.5897 |
| 0.05 | 7568 | 0.5899 |

Centre slip moves 0.07 mm over the eightfold sweep. Per-column max absolute
difference against `xi = 0.2` across all 9 stations: slip_2 0.74 %, V2 0.33 %,
tau_2 0.60 %, p 0.42 %, state 1.70 %. All under 2 %. `xi = 0.05` costs 42.8 %
more steps for no resolvable gain, so `xi = 0.2` is the right setting.

### Domain — **not** converged, and this gold sits at the worst corner

| x, z | y | grading | centre (mm) | edge (mm) | edge/centre |
|---|---|---|---|---|---|
| **±500** | **±500** | 1.15 | **37.31** | **21.99** | **0.5895** ← this gold |
| ±1000 | ±500  | 1.15 | 36.94 | 19.58 | 0.5302 |
| ±1500 | ±500  | 1.15 | 36.90 | 19.52 | 0.5291 |
| ±500  | ±2000 | 1.15 | 38.33 | 29.66 | 0.7739 |
| ±500  | ±2000 | 1.0  | 38.32 | 29.44 | 0.7682 |

In x and z the ratio change shrinks 59-fold across the two refinement steps
(0.059 then 0.001), which is **consistent with convergence but not established
from two points** — the ±3000 run that would confirm it is still running. In y
it is plainly not converged: moving the rigid wall from ±500 to ±2000 shifts the
ratio 31 % with no sign of saturation, and only two y points have been tested.
Mesh grading (`enlarging_ratio` 1.15 against 1.0) changes it by 0.7 %.

### Dip-direction slip is real off the symmetry axes

Five of the nine stations lie on `strk = 0` or `dp = 0`, where antiplane
symmetry forces dip slip to zero; there `slip_3` is 1e-19 or smaller and V3 sits
pinned at the `V_zero = 1e-20` floor, swinging on round-off. **The other four do
not.** At `strk = ±200` *and* `dp = ±200`, `slip_3` is a stable 2.89e-4 m —
1.89 % of the co-located `slip_2`, converged to under 1 % across the whole `xi`
sweep — with V3 at −10.68 and an antisymmetric sign pattern.

That is physical: a slipping patch in 3-D produces shear-stress changes with
both components, and only the symmetry axes force the dip component to zero.
Section 4.1 asks for V3, and at these four stations the submitted column carries
signal.

Consequence for anyone comparing runs here: a relative-difference check
normalised per column reports meaningless percentages for the five on-axis
stations, so **any relative check needs an absolute floor** — and must not then
conclude the whole column is noise.

### What the study establishes, and what it does not

**Establishes:** `xi = 0.2` and `dx = 50 m` are adequate; the y extent dominates
the domain error; grading is not a factor.

**Does not establish:** a converged x, z extent (pending ±3000); a converged y
extent (two points, both far from flat); or that this gold's ±500 configuration
is right — it is the corner the sweep argues against. Its agreement with Kim
(37.31 vs ~38 mm centre, 21.99 vs ~21 edge) is partly luck; the better-resolved
domains give 36.9 mm and edge/centre 0.53.

## Regenerating

Run `case_input/test.bp8.qdc` with `xi = 0.2`, `nstep = nt_out = 8000`, which
exits on `fluid_tend` at 30.0 days. Anything less than ~5301 steps is a
truncated run: the compset default `nstep = 200` reaches only 1.14 days, and
comparing that against this oracle produces a failure that is not real.
