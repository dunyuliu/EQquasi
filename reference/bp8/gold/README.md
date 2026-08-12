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

## Regenerating

Run `case_input/test.bp8.qdc` with `xi = 0.2`, `nstep = nt_out = 8000`, which
exits on `fluid_tend` at 30.0 days. Anything less than ~5301 steps is a
truncated run: the compset default `nstep = 200` reaches only 1.14 days, and
comparing that against this oracle produces a failure that is not real.
