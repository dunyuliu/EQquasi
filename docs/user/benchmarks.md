# Benchmarks

EQquasi is exercised on problems from the SCEC Sequences of Earthquakes and
Aseismic Slip (SEAS) benchmark project
(https://strike.scec.org/cvws/seas/benchmark_descriptions.html) and on a
few problems of its own. Results from v1.2.1 on SEAS BP5 are published in
Jiang et al. (2022, JGR); see [Citing](citing.md).

## How results are checked

A set of frozen results lives under `data/<compset name>/`, each named
for the compset that produced it. The test suite reruns those compsets and
compares every output file against the reference.

These references are **regression locks, not validations**. They detect an
unintended change in the results; they do not by themselves establish that
a benchmark is reproduced correctly. The comparison is relative, never
bit-for-bit: two runs of the same case on the same host differ by about
1e-14 from the order of MPI reductions, so each value is judged against the
size of the quantity it belongs to.

```
python3 -m pytest tests/              # source and format checks, no runs, about 1-2 minutes
python3 -m pytest tests/ -m e2e_fast  # builds and runs the small cases, about 20 minutes
python3 -m pytest tests/ -m e2e       # adds full first cycles, about 75 minutes and longer
```

The e2e tiers need EQquasi installed for the current source
(`bin/eqquasi-<version>`), with `EQQUASIROOT` set and `scripts/` on `PATH`.

## SEAS benchmark cases

| benchmark | compset | resolution | what it is |
|---|---|---|---|
| [BP5](https://strike.scec.org/cvws/seas/benchmark_descriptions.html) | `bp5.qdc.2000` | 2000 m | vertical strike-slip fault, rate-and-state friction with the aging law, 3D |
| [BP7](https://strike.scec.org/cvws/seas/download/SEAS_BP7_May2023rev.pdf) | `bp7.qdc.a.10` | 10 m | small velocity-weakening patch with a prescribed nucleation perturbation |
| [BP8](https://strike.scec.org/cvws/seas/benchmark_descriptions.html) | `bp8.qdc.gs.10` | 10 m | fluid injection with a Gaussian source and along-fault pore-pressure diffusion |

Their test versions -- `test.bp5.qdc.2000`, `test.bp7.qdc.a.10` (at 25 m)
and `test.bp8.qdc.gs.10` (at 50 m) -- are the ones the test suite runs.

**BP7 nucleates off-centre on purpose.** The benchmark prescribes a shear
traction perturbation centred at (-50 m, -50 m), ramped in over 1 s. A BP7
rupture that nucleates at the origin would be wrong. The perturbation is
built into the solver (`bp = 7`), not set in the compset.

**BP5's nucleation patch starts out of steady state on purpose.** The patch
is given a high slip rate and a state variable set for the background creep
rate, and that mismatch is what makes it accelerate. Shear stress is set
for the node's own slip rate.

## Compsets

Every compset, with whether a reference result checks it:

| compset | what it is | checked by |
|---|---|---|
| `bp5.qdc.2000` | SEAS BP5 | full first cycle, via `test.bp5.qdc.2000` |
| `bp5.qdc.kink.2000` | BP5 friction, fault with a 10 degree bend | full first cycle |
| `bp7.qdc.a.10` | SEAS BP7 | via `test.bp7.qdc.a.10` |
| `bp8.qdc.gs.10` | SEAS BP8, Gaussian source | via `test.bp8.qdc.gs.10` |
| `bp1002.qdc.2500` | BP5 friction, two faults with a 5 km step-over | full first cycle |
| `bp1002.qdc.caps.2500` | `bp1002.qdc.2500` with normal-stress caps | not checked |
| `bp1002.qdc.caps.taper.2500` | with caps and a 5 km velocity-weakening taper at the inner tips | not checked |
| `bp1002.qdc.zone.2500` | step-over with unequal segments and a fixed velocity-weakening patch | not checked |
| `liu2020.qdc.kink.600` | bent fault of Liu, Duan and Luo (2020) at 600 m | reference result, no automatic run |
| `liu2020.qdc.kink.300` | the same at the paper's 300 m | not checked |
| `das.qdc.10` | older quasi-dynamic case at 10 m | not checked |
| `enechelon.qdc.gap4.1000` | two en echelon faults, 4 km step-over | not checked |
| `bp1001.fdc.250` | fully dynamic mode, planar fault | not checked |
| `bp1001.fdc.rough.250` | fully dynamic mode, rough fault | not checked |
| `bp1001.qdc.rough.250` | quasi-dynamic, rough fault | not checked |
| `liu2020.fdc.planar.300` | fully dynamic mode, planar fault | not checked |
| `liu2020.fdc.rough.250` | fully dynamic mode, rough fault | not checked |

The `fdc` compsets set `mode = 2`, in which a cycle stops as soon as a
rupture begins, so a dynamic rupture code can take over; `qdc` compsets
(`mode = 1`) stop after the rupture has ended.

"Not checked" means no reference result exists, so nothing tests the
numbers; it does not mean the compset is broken. `liu2020.qdc.kink.300` and
`.600` do not reproduce the paper: the paper's ruptures come from the
dynamic code EQdyna, and this is EQquasi alone. At 600 m the cohesive zone
is under-resolved.

Test versions, small and fast and not for science:

| test compset | what it is |
|---|---|
| `test.bp5.qdc.2000` | BP5, 101 steps; a full cycle in the long tier |
| `test.bp5.qdc.dip90.2000` | BP5 friction, planar control for the kink case, 101 steps |
| `test.bp7.qdc.a.10` | BP7 at 25 m, 101 steps; a full cycle in the long tier |
| `test.bp8.qdc.gs.10` | BP8 at 50 m, run to 30 days |
| `test.stepover.qdc.1000` | two faults, releasing step-over, 101 steps |
| `test.stepover.qdc.con.1000` | two faults, restraining step-over, 101 steps |
| `test.bp1002.qdc.zone.2500` | test version of the zoning experiment; no reference yet |

## Submitting BP8 results

BP8 results go to the CRESCENT DET platform as a zip. `checkBP8Submission`
checks a result directory against the benchmark description and writes the
zip:

```
checkBP8Submission result/cycle0 --zip auto
```

`--zip auto` writes the zip next to the result directory and names it
`<modeler>_eqquasi-<version>-<dx>m.zip`, with the version (dots turned into
dashes) and cell size read from `runInfo.json`; `--modeler` sets the first
part. It refuses to write a zip while any check fails. Section 4.3 wants profile nodes every 10 m exactly, so a run at a
coarser `dx` must first be resampled onto that grid (the interpolation is
recorded in each file header):

```
resampleBP8Profiles.py result/cycle0 submission
checkBP8Submission submission --zip auto
```

Two settings matter for a comparable BP8 run. `nstep` must be large enough
to reach `fluid_tend` (30 days): at `xi = 0.2` that takes about 5300 steps,
and a smaller cap stops the run early. One MPI rank with
`OMP_NUM_THREADS=1` is the sensible choice; the problem is small and extra
ranks barely help.

## Further reading

Each reference directory has a README with what that result does and does
not establish. The BP8 one (`data/test.bp8.qdc.gs.10/README.md`)
covers the pore-pressure solver, the domain-size and time-step studies, and
the initial condition.
