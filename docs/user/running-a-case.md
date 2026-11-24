# Running a Case

## Create a case

```
create.newcase scratch/mycase bp5.qdc.2000
cd scratch/mycase
```

The second argument is a compset: a predefined case, one directory under
`compset/` holding a `user_defined_params.py`. The full list, with what each
one is for, is on the [Benchmarks](benchmarks.md) page. Compset names follow
`<benchmark>.<mode>[.<variant>].<dx_m>[.<description>]`, where `mode` is
`qdc` (quasi-dynamic) or `fdc` (fully dynamic) and `dx_m` is the on-fault
element size in metres. A `test.` prefix marks a small, short version used
by the test suite, not meant for science.

Keep cases under `scratch/` in the repository (it is gitignored) or anywhere
outside the repository. `create.newcase` deletes and replaces an existing
directory of the same name.

## What a case holds

```text
<case>/
    user_defined_params.py   the only file you edit
    case.setup  case.submit  run.sh    what you run
    input/                   everything a cycle reads; frozen once a cycle has run
    result/                  cycle0/ ... cycleN/, one per completed cycle
    scratch/                 the cycle currently running, and nothing else
    tool/                    modules the case imports, not for editing
    log/                     setup and per-cycle run logs
```

## Configure

Edit `user_defined_params.py`. Every setting is an attribute of `par`; see
[Parameters](parameters.md) for all of them. The ones most often changed:

* `istart`, `iend` -- the earthquake cycles to run, numbered from 1.
* `HPC_ncpu` -- the number of MPI ranks `run.sh` launches.
* `nstep` -- the most time steps a cycle may take before it is stopped.
* `nt_out` -- how often (in steps) fault and volume snapshots are written.
* `dx` -- the on-fault element size, in metres.

### Time integrator and final time

`time_integrator` selects how a cycle is advanced in time.

* `1` (the default) is the loop every reference result was made with. Each
  step takes an explicit (Euler) predictor of the fault displacement, a
  trapezoidal corrector, a Newton solve for the slip rate after each of the
  two, and then one exponential update of the state at the corrected slip
  rate. The step size is `xi*Dc/Vmax` (per node when `xi_gamma > 0`), rounded
  down to a whole number of base steps `dt`, never below `dt`, and capped by
  `dtmax`.
* `2` is an adaptive Runge-Kutta integrator, the Dormand-Prince 5(4) pair
  (Dormand and Prince, 1980), experimental. Its unknowns are the fault
  displacement and the state (ln(V0*theta/Dc) for `friclaw = 3`, psi for
  `friclaw = 4`); the slip rate comes from the same Newton solve at every
  stage. An attempted step costs 6 slip-to-traction products (value 1: 2 per
  step) and is accepted or repeated with a smaller size from its own error
  estimate. `xi` and `xi_gamma` are not applied and the step may fall below
  `dt`; `dtmax` still caps it. The first step size comes from the initial
  rates (Hairer, Norsett and Wanner, 1993). The tolerances are fixed in the
  source (relative 1e-6; absolute 1e-9 m on displacement and 1e-6 on the
  state); they are not calibrated yet and are not parameters.

Value 2 runs with `solver = 4`, `bp = 5` or `7`, `friclaw = 3` or `4` and
`fric_pc_L = 0` only, and only the first cycle: it cannot continue from a
previous cycle, and, as with value 1, a run that was killed starts again from
step 0. Any other combination stops with code 511; nothing falls back to
value 1. BP8 is refused because its pore pressure is advanced once per step.

With value 2 the cycle log holds one line per accepted step,
`RK45 step= ... t= ... h= ... dtev= ... err= ... Vmax= ... rej_step= ... rej_total= ...`
(`h` the step taken, `dtev` the `xi*Dc/Vmax` value for comparison, `err` the
error norm, at most 1), and one `RK45 summary` line with the accepted and
rejected steps and the number of products. The `pma` value printed on screen
is not meaningful with value 2.

`tmax` (seconds, default 0 = off) stops a run at that simulated time; the last
step is cut to end on it. `nstep` and the slip-rate exit still apply if they
come first. It is honoured by both integrators with `solver = 3` or `4`
(`solver = 1` or `2` stops with code 511), so two runs can be compared at one
time instead of one step count.

## Set up

```
./case.setup
```

`case.setup` writes the solver inputs (`model.txt`, `stations.txt` and
`on_fault_vars_input.nc`) into `input/`, and writes `run.sh`. A compset with
a rough or bent fault also has `bFault_Rough_Geometry.txt` there. It also prints checks worth reading before a long run:

* the estimated mesh size and resource need;
* the cohesive-zone width against the cell size -- below about one cell the adaptive time step can oscillate instead of resolving a rupture;
* whether the largest allowed time step can overload the fault in one step;
* whether normal-stress caps are on, with a warning when a bent, rough or multi-fault geometry runs without them.

Once a case has run, `case.setup` refuses to run again: regenerating
`input/` would give later cycles different inputs from earlier ones. Start
a new case for changed parameters. If you are certain the regenerated
inputs are identical, `./case.setup --force` overrides the refusal.

`case.setup` is copied into the case when it is created, so a case keeps the
version it was made with. To update an existing case to the current one,
run `cp $EQQUASIROOT/scripts/case.setup .` in the case.

### Compsets with a geometry generator

Compsets with a bent fault (`bp5.qdc.kink.2000`, `liu2020.qdc.kink.300`,
`liu2020.qdc.kink.600`) ship a script that writes the fault surface. Run it
after `create.newcase` and before `case.setup`, from the case directory:

```
python3 input/generateKinkGeometry.py input
./case.setup
```

## Run

```
bash run.sh
```

`run.sh` loops over cycles `istart` to `iend`. Each cycle runs inside
`scratch/`, its screen output is copied to `log/cycle<N>.log`, and only when
it succeeds is `scratch/` moved to `result/cycle<N>/`. The next cycle
restarts from that cycle's `fault.r.nc` and `disp.r.nc`. Cycle `i` writes
`result/cycle<i-1>/`, so `istart = 1` produces `result/cycle0/`.

A cycle ends when the fault has ruptured and the peak slip rate has stayed
below `exit_slip_rate` for 10^5 s of simulated time, or when it reaches
`nstep` steps. BP8 cases instead stop at the simulated time `fluid_tend`.

If a cycle fails, `run.sh` stops, prints the exit status and the last 20
lines of the log, and leaves `scratch/` as it was for inspection. No
`result/cycle<N>/` is created for a failed cycle, so a half-written cycle
cannot be mistaken for a finished one. To retry, remove the failed attempt
and run again:

```
rm -rf scratch
bash run.sh
```

To continue a finished case with more cycles, set `istart` to the next
cycle and raise `iend` in `user_defined_params.py`, then run
`./case.setup --force` and `bash run.sh`. Those two values change only
`run.sh`, so the regenerated inputs are identical to the ones already used.

### MPI launcher and ranks

`run.sh` starts `HPC_ncpu` MPI ranks with the launcher recorded when
EQquasi was installed. Three environment variables override it for one run
without editing anything:

```
MPIRUN=/path/to/mpirun bash run.sh
MPIRUN_ARGS="--bind-to none" bash run.sh
OMP_NUM_THREADS=1 bash run.sh
```

* `MPIRUN` -- the launcher. It must come from the same MPI the binary was built against.
* `MPIRUN_ARGS` -- extra launcher flags. The installer records a default: `--bind-to none` for Open MPI, plus shared-memory transport flags on conda builds, and none for `ibrun`.
* `OMP_NUM_THREADS` -- threads per rank; `run.sh` defaults it to 1.

With Open MPI's default core binding, two runs started on the same host are
both pinned to the same first cores and each runs at about half speed;
`--bind-to none` lets the operating system spread them. Extra ranks help
less than you might expect, because the direct solve does not distribute
well (see [Performance](performance.md)).

### On a cluster

On TACC Lonestar6 the installer records `ibrun` as the launcher. Submit the
case with

```
./case.submit
```

which runs `sbatch batch.hpc`. `batch.hpc` holds the SLURM header built from
the `par.HPC_*` parameters (nodes, ranks, queue, wall time, account, email),
loads the modules, and runs `bash run.sh`, so a batch job runs exactly what
an interactive `bash run.sh` would. You can also run `bash run.sh` yourself
from inside an allocation. `batch.cycle.eqquasi.hpc` belongs to the coupled
EQquasi/EQdyna cycle workflow and still uses the older flat case layout.

## Plot

Post-processing tools live in `scripts/` and are on your `PATH` once
`install.eqquasi.sh` has been sourced. Run them from inside the case. Each
takes cycle directories as arguments, or none to process every cycle:

| tool | writes | scope |
|---|---|---|
| `plotRuptureTime.py` | `rupture_time.png` | per cycle |
| `plotPeakSliprateTime.py` | `peak_slip_rate_vs_time.png`, `peak_slip_rate_per_fault_vs_time.png` | all cycles |
| `plotOnFaultVars` | `fault.NNNNN.nc.png` per snapshot, `on_fault_vars.gif` | per cycle |
| `plotAccumulated` | `accumulatedSlip.horizontal.png` (or `.vertical`) | all cycles, stacked |
| `plotStations.py` | `on_fault_stations.png`, `off_fault_stations.png` | per cycle |
| `plotMagnitudeTime.py` | `mag_vs_time.png`: moment magnitude of each event against time | all cycles |

```
plotRuptureTime.py result/cycle2
plotPeakSliprateTime.py result/cycle0 result/cycle1 result/cycle2
plotAccumulated --depth-km -10 result/cycle0 result/cycle2
plotStations.py result/cycle2
plotMagnitudeTime.py
```

Each tool's `-h` lists its options. `-o DIR` sends the figures elsewhere.

**Multi-fault runs.** `plotAccumulated` and `plotRuptureTime.py` draw one
fault per figure, chosen with `--fault` (counted from 0); plot every fault, and pass the same
`--ylim` to each so the figures are comparable. `plotOnFaultVars` loops over
the faults itself and tags its figures `.f0`, `.f1`, and so on.
`plotPeakSliprateTime.py` draws the whole model in one figure and each fault
in a second.

**Reading accumulated slip.** The profile crosses the velocity-weakening
zone, the velocity-strengthening margins and the creeping region outside
rate-and-state control, so the largest value on the plot can be imposed
creep rather than earthquake slip. Restrict to the velocity-weakening zone
before quoting a coseismic slip. Each cycle's curve starts from the sum of
the cycles before it.
