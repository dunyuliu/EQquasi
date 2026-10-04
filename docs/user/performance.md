# Performance

Every figure on this page is taken from a measurement recorded in the
repository: the `runInfo.json` of each reference result, and timing studies
in the BP8 reference notes. They were measured on shared workstations, on
the versions named, and are a guide to relative cost rather than a promise
for your hardware. Each run's own `runInfo.json` records its timings, so
you can compare your runs the same way.

## Reference runs

From `runInfo.json` in each `reference/` directory. "Time loop" excludes
mesh generation, assembly and the first factorization; `OMP_NUM_THREADS=1`
throughout.

| case | version | CPU | ranks | elements | steps | time loop (s) | s/step |
|---|---|---|---|---|---|---|---|
| `test.bp5.qdc.2000`, full cycle | 1.7.2 | AMD EPYC 7532 | 4 | 68,400 | 4483 | 3087 | 0.69 |
| `test.bp5.qdc.2000`, 101 steps | 1.7.0 | AMD EPYC 7532 | 2 | 68,400 | 101 | 85 | 0.84 |
| `test.bp5.qdc.dip90.2000` | 1.16.0 | AMD EPYC 7543 | 2 | 68,400 | 101 | 78 | 0.77 |
| `bp5.qdc.kink.2000`, full cycle | 1.14.0 | AMD EPYC 7543 | 3 | 68,400 | 4876 | 3177 | 0.65 |
| `test.bp7.qdc.a.10`, full cycle | 1.7.0 | AMD EPYC 7532 | 2 | 64,000 | 1358 | 1013 | 0.75 |
| `test.bp7.qdc.a.10`, 101 steps | 1.7.0 | AMD EPYC 7532 | 2 | 64,000 | 101 | 80 | 0.79 |
| `test.bp8.qdc.gs.10`, 30 days | 1.6.0 | AMD EPYC 7532 | 1 | 8,000 | 5301 | 796 | 0.15 |
| `bp1002.qdc.2500`, full cycle | 1.11.0 | AMD EPYC 7532 | 3 | 41,472 | 3821 | 2615 | 0.68 |
| `test.stepover.qdc.1000` | 1.16.0 | AMD EPYC 7543 | 2 | 52,500 | 101 | 80 | 0.79 |
| `test.stepover.qdc.con.1000` | 1.19.0 | AMD EPYC 7F72 | 2 | 52,500 | 101 | 79 | 0.78 |
| `liu2020.qdc.kink.600`, 10000 steps | 1.14.0 | AMD EPYC 7543 | 3 | 210,000 | 10000 | 20,254 | 2.03 |

For a whole BP5 first cycle at 2000 m, plan on about 27 minutes on TACC
Lonestar6 or about 51 minutes on 4 ranks of a shared workstation.

## Parallel scaling

Measured on the 8000-element `test.bp8.qdc.gs.10` case, 200 steps,
`OMP_NUM_THREADS=1`, on an otherwise idle machine:

| MPI ranks | wall (s) | s/step | speedup |
|---|---|---|---|
| 1 | 29.4 | 0.147 | 1.00x |
| 2 | 25.6 | 0.128 | 1.15x |
| 4 | 24.1 | 0.121 | 1.22x |
| 8 | 23.8 | 0.119 | 1.24x |

Eight times the ranks buys 24%. At this problem size the direct
factorization is dominated by work that does not distribute, so add ranks
mainly when the mesh is large enough that memory, not time, is the limit.
Always set `OMP_NUM_THREADS` (`run.sh` sets it to 1 unless you do): the
build uses OpenMP, and unpinned threads from several ranks oversubscribe a
node.

## Scaling with mesh size

Two `test.bp8.qdc.gs.10` runs differing only in resolution:

| dx | elements | s/step | peak memory (RSS) |
|---|---|---|---|
| 50 m | 12,800 | 0.196 | 234 MB |
| 25 m | 64,000 | 0.810 | 1.23 GB |

Five times the elements cost about 4.1 times the time per step and 5.2 times
the memory, close to linear over this range. Beyond it, treat any
extrapolation as an order of magnitude only: the cost of a direct sparse
factorization eventually grows faster than the number of unknowns.

## Measuring

Timings on a shared machine depend on what else is running. One `dx = 25 m`
timing on the same hardware came out ten times slower because the machine
was under a load average of about 33. Check `uptime` before timing a run.

## Hardware

The workstation figures above come from AMD EPYC 7532 (64 logical cores),
EPYC 7543 and EPYC 7F72 Linux workstations, as recorded in each run's
`runInfo.json`. No GPU path exists. Systematic scaling on Lonestar6 or other
clusters has not been measured.
