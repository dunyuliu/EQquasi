# BP5 gold, cycle 0

The first full earthquake cycle of `case_input/test.bp5.qdc` — interseismic
loading, nucleation, rupture, and postseismic decay — kept alongside the
step-101 snapshot in the parent directory rather than replacing it.

## Why both

`../fault.00101.nc` stops at step 101, before nucleation. It is a tight lock
(the four gates compare it at max abs diff **0.0**) but it only constrains the
initial transient: a change that altered how the fault ruptures would not move
it at all.

This directory constrains the part that matters. 310 of 1891 fault nodes
rupture, over 0.065–117 s, peaking at 1.16 m/s.

## Configuration

| | |
|---|---|
| compset | `test.bp5.qdc`, `dx = 2000 m` (BP5's production resolution) |
| `nstep` | 100000 — high enough that the cycle exits on its own criterion, not a step cap |
| `nt_out` | 1000 |
| exit | `exit_slip_rate = 1e-3`, then the hang-time window below |
| run | 4483 steps, 1.166 days simulated |
| cost | 3087 s on 4 MPI ranks, 0.689 s/step; factorization 4.4 s |
| mesh | 75640 nodes, 68400 elements, 204228 equations, 1891 fault nodes |

## Files

| | |
|---|---|
| `fault.0*.nc` | on-fault snapshots every `nt_out`, plus the final step |
| `fault.r.nc` | restart state at exit |
| `disp.0*.nc`, `disp.r.nc` | volumetric displacement field at the same cadence |
| `fltst_*.txt` | 9 on-fault station time series |
| `srfst_*.txt` | 9 off-fault station time series |
| `global.dat` | 4483 × 7 |
| `cplot_EQquasi.txt` | fault-plane final state, 1891 × 16; column 16 is `fnft`, rupture time |
| `cplot_ruptarea_trac_slip.txt` | rupture area, traction, slip |
| `tdyna.txt` | rupture start and end times |
| `runInfo.json`, `summary.json` | provenance and the scalars worth asserting |

17 MB. The `disp.*.nc` files are most of it — the full field over 75640 nodes.
If this pattern is repeated for BP5-dip90 and BP7, consider keeping only the
on-fault snapshots, which is what every comparison actually reads.

## The one thing to know before trusting it

**The cycle length is set by a magic number in the source.** `exitCriteria`
(`src/solveTimeLoopMUMPS.f90`) does not exit when max slip rate falls below
`slipr_thres`; it records `tdynaend` and then waits for the low-slip-rate state
to persist for `hang_time = 1e5` s of simulated time, hardcoded with the comment
"the 100 seconds threshold is subjective to change". Change it and this oracle
moves, through no change in physics. It belongs in `user_defined_params.py`.

Like every reference here, this is a regression lock, not a validation: it
detects unintended change, it does not establish that BP5 is reproduced
correctly.
