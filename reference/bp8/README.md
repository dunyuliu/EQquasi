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

### BP8 pore fluid diffusion

The along-fault pore pressure solver used by `bp == 8` is verified against the
closed-form solution in eq. (22) of the SEAS BP8 benchmark description, which
gives `p(x2, x3, t)` for a Gaussian injection source in an unbounded domain.

Relative error in pore pressure at the injection point, `BP8-QD-GS`:

| t (s)   | `dx` = 50 m | `dx` = 25 m |
|---------|-------------|-------------|
| 25 000  | 5.11 %      | 1.51 %      |
| 50 000  | 5.75 %      | 1.50 %      |
| 75 000  | 5.74 %      | 1.44 %      |
| 100 000 | 5.60 %      | 1.37 %      |

Halving `dx` reduces the error by a factor of about 3.9, i.e. **second-order
convergence**, as expected for the FTCS Laplacian combined with a discretely
normalized Gaussian source. Extrapolating to the production resolution
`dx = 10 m` gives roughly 0.2 %.

The remaining error is source under-resolution, not a modelling error: the
Gaussian has a characteristic size `L_gauss = 50 m`, so at `dx = 50 m` it is
resolved by a single cell and at `dx = 25 m` by two.

Two caveats on the comparison itself:

  - The analytic solution assumes an unbounded fault, whereas the solver imposes
    zero fluid flux on the boundary of the frictional domain at +/-400 m. The
    two only agree while the diffusion length `sqrt(4*alpha*t)` stays well
    inside that boundary; it reaches 400 m at around `t = 8e5 s`, i.e. roughly
    9 days into the 30-day benchmark.
  - Eq. (20) of the benchmark description writes the Gaussian source as
    `(1/L^2) exp(-r^2 / 2L^2)`, which integrates to `2*pi` rather than 1, and is
    therefore inconsistent with the closed form in eq. (22) by that factor. The
    implementation follows eq. (22)/(25) and normalizes the source weights
    discretely so that they sum to exactly 1.

### Elastic domain size matters, and this section used to say otherwise

An earlier version of this text concluded, from an eightfold increase in
elements at fixed extent, that "the elastic domain size does not matter here".
That was wrong, and it was wrong in a way worth recording: the study varied
*resolution* while holding the *extent* fixed, so it could not see the effect it
claimed to rule out.

The domain sweep in this file's "Domain truncation" section shows edge/centre
slip moving from 0.589 to 0.529 as x and z go from +-500 m to +-1500 m, and to
0.774 when the y wall moves to +-2000 m. Extent matters more than anything else
measured on this benchmark.

### The time-step factor xi can be loosened

`xi` is the Lapusta et al. (2009) safety factor in `dtev = xi * Dc / V_max`
(`par.xi`). BP5 uses 0.015; BP7 and BP8 use 0.05. At `dx = 50 m`, 4000 steps
each:

| xi   | wall   | simulated time reached | peak log10 V | peak moment rate | slip     |
|------|--------|------------------------|--------------|------------------|----------|
| 0.05 | 9m22s  |  6.03 d                | -6.377       | 1.029e9          | 42.76 mm |
| 0.1  | 9m22s  | 17.00 d                | -6.376       | 1.032e9          | 42.74 mm |
| 0.2  | 9m27s  | 22.33 d                | -6.374       | 1.032e9          | 42.74 mm |

Same wall clock, 3.7x more simulated time, and the answers agree to 0.003 log
units in peak slip rate, 0.3 % in moment rate and 0.05 % in slip. For BP8 at
this resolution `xi = 0.2` is essentially free, cutting a 30-day run from about
20 000 steps to about 5 500.

Note this was measured at `dx = 50 m`. A finer mesh resolves faster transients,
so re-check before assuming it carries to `dx = 10 m`.

### Submitting to the CRESCENT DET platform

```
checkBP8Submission <result_dir> --zip auto --both-station-names
```

validates a result directory against the benchmark description and writes the
upload zip. Naming convention:

```
<modeler>_<code>-<version>-<resolution>m.zip
e.g. dliu_eqquasi-1-4-3-10m.zip
```

Version and resolution are read from `runInfo.json`, not typed, so the dataset
label on the platform always matches the binary and mesh that produced the data.
Version dots become dashes so everything after the underscore is one
dash-separated token, which is what the platform displays.

All BP8 output files use the `.dat` extension, **including the station time
series** (`fltst_strk+000dp+000.dat`). The platform routes purely on filename
and processes `.dat` for every file type; a station file written as `.txt` is
silently ignored, so the section 4.3 profiles render and the section 4.1 time
series do not, with no error reported anywhere. Section 4.1 of the description
lists station names bare, without any extension, which is what led us to `.txt`
in the first place -- the platform's own `processed_files` listing is the
authority.

BP5 and BP7 keep `.txt`, which is their existing convention.

### Full 30-day run and discrete mass conservation

`test.bp8.qdc` run to the benchmark's full `t_f = 30 days` (5185 steps at
`dx = 50 m`, 17 min on 2 ranks) exits on the time criterion as intended and
gives an independent check that does not rely on the Green's function at all.

After injection stops at `t_off = 100 h`, the zero-flux condition on the
boundary of `Omega_f` traps every injected drop inside the frictional domain,
so the pore pressure must relax towards a uniform value fixed only by the total
volume injected, `V = Q0 * t_off = 1080 m^3`:

```
p_uniform = V / (beta * phi * A * L_fwid)
```

| quantity                                | value       |
|-----------------------------------------|-------------|
| predicted, discrete area (17x17 nodes)  | 1.4948 MPa  |
| **code, centre station at t = 30 d**    | **1.5013 MPa** |
| error                                   | **+0.43 %** |

This exercises the diffusion operator, the source normalization and the
boundary condition together.

**A discretization detail this exposes.** The agreement holds against the
*discrete* area, `(17 * 50)^2 = 722 500 m^2`, not the continuum area of
`Omega_f`, `800^2 = 640 000 m^2`. Pore pressure lives at nodes and every node
carries a full `dx^2` cell, including those on the boundary, so the effective
storage area overshoots by one cell width on each side. At `dx = 50 m` that is
a 13 % area overestimate and therefore an 11 % *underestimate* of the
long-time uniform pressure. The error is first order in `dx` and shrinks to
about 2.5 % at the production `dx = 10 m`, but it is a real bias in late-time,
boundary-influenced results and is not removed by the second-order convergence
reported above, which is measured while the pressure front is still far from
the boundary.

Physical behaviour over the full run, at the centre station:

From `work/bp8.prod10`, dx = 10 m, `xi` = 0.2, literal specification initial
condition (see "The initial condition is over-determined" below):

| quantity                        | value                              |
|---------------------------------|------------------------------------|
| peak pore pressure              | 13.051 MPa at `t_off` = 4.17 d     |
| minimum effective normal stress | 11.949 MPa (52 % reduction)        |
| peak slip rate                  | 3.2e-7 m/s (log10 -6.489) at 0.80 d |
| secondary slip-rate peak        | log10 -6.804 at 4.37 d             |
| total slip at 30 d              | 45.04 mm                           |

Effective normal stress never reaches zero.

An earlier revision of this section reported a peak at 8.6 d and read it as "the
expected signature of a slip front that keeps propagating once injection ends".
That reading was wrong. The dominant peak is at **0.80 d**, well *before*
shutoff, and it is an artefact of the over-determined initial condition: the
patch starts 1.67 MPa overstressed and sheds that stress as an early aseismic
transient. The physical response to shutoff is the *secondary* peak at 4.37 d.
Removing the overstress removes the early peak entirely and leaves a maximum of
about log10 -6.80, which is close to the amplitude of the secondary peak here.

The friction solve was checked independently of the fluid: with
`tau0 = 14.6 MPa` on `sigma_bar = 25 MPa` and `theta = Dc/V_init`, the
regularized rate-and-state law gives `V = 6.54e-11 m/s`, and the code reports
`6.55e-11 m/s` at the first output step.

That check was narrower than it looked: it confirms the solver inverts the
friction law correctly, but not that `6.54e-11 m/s` is the value the benchmark
wants. It is not what eq. (27) asks for, and that leads to the following.

### Running BP8-QD-GS: the recipe

```
export EQQUASIROOT=$PWD PATH=$PWD/bin:$PWD/scripts:$PATH OMP_NUM_THREADS=1
EQQUASIROOT=$PWD MACHINE=utig make -C src && mv src/eqquasi bin/

create.newcase work/mycase bp8.qdc.gs.10     # dx = 10 m production
#              work/mycase test.bp8.qdc      # dx = 50 m test
cd work/mycase && ./case.setup && echo 1 > currentcycle.txt
setsid nohup mpirun -np 1 $EQQUASIROOT/bin/eqquasi > r.log 2> t.log < /dev/null &
```

Then package for the CRESCENT DET platform:

```
python3 scripts/resampleBP8Profiles.py work/mycase work/mysub   # only needed if dx > 10 m
python3 scripts/checkBP8Submission work/mysub --zip auto
```

Points that are easy to get wrong, each of which cost real time here:

  - **`np = 1`, `OMP_NUM_THREADS = 1`.** The problem is small; extra ranks and
    threads do not pay for themselves.
  - **`nstep` must be large enough to reach `fluid_tend`.** The run exits on
    `t_f` = 30 days, but only if `nstep` allows it -- at `xi` = 0.2 and
    `dtmax` = 500 s that needs about 5300 steps. Capping `nstep` at 4000 stops
    the run at ~23 d and quietly produces a short, non-comparable result.
  - **Only dx = 10 m is directly submittable.** Section 4.3 wants profile nodes
    at exactly 10 m over [-400, 400]; a coarser run must go through
    `resampleBP8Profiles.py`, which interpolates onto that grid and says so in
    the headers.
  - **Never `pkill -f <pattern>`** to stop a run. The pattern matches the shell
    issuing it, so it kills itself; that destroyed two long runs here. Kill by
    PID, found via `readlink /proc/<pid>/cwd`.

### Elastic domain truncation

Section 6 of the description warns that truncating the whole-space "will most
likely change results at least quantitatively, or even qualitatively" and asks
that boundaries be extended until results stop moving. That warning is
justified: `meshgen.f90` constrains only the two fault-normal faces, so the `x`
and `z` faces are traction-free, and at the default `+-500 m` they sit just
100 m outside `Omega_f`.

Measured at dx = 50 m, 30 days, both initial-condition readings:

| run | box half-width | slip `(0,0)` | slip `(-200,0)` | edge/centre |
|-----|----------------|--------------|-----------------|-------------|
| reading A | 500 m  | 42.78 mm | 28.79 mm | 0.673 |
| reading A | 2000 m | 41.63 mm | 25.47 mm | 0.612 |
| reading C | 500 m  | 37.31 mm | 21.99 mm | 0.589 |
| reading C | 2000 m | 36.90 mm | 19.52 mm | 0.529 |
| `taehoKim_ref` | (BEM, none) | ~38 mm | ~21 mm | 0.553 |

Two things follow. Enlarging the box barely moves the amplitude at the injection
point (-1 to -3 %) but moves the slip at 200 m substantially (-11 %), so
truncation controls the *shape* of the slip profile rather than its size. And
the reference's 0.553 falls **between** our 500 m and 2000 m results, so we now
bracket it rather than sitting on one side; the 500 m box happens to match the
two amplitudes best (-2 % and +5 %) while the 2000 m box matches the shape
better. Neither is converged, and that is worth stating rather than presenting
the closer one as agreement.

An earlier sweep in this project concluded the results were domain-independent.
It compared only centre-station slip and peak `Vmax` -- the two quantities least
sensitive to truncation -- and so measured the wrong thing. A cause is not
excluded until the quantity it would actually affect has been measured.

### The initial condition is over-determined

Equation numbers below follow the **2026-08-06 revision** of the description,
which dropped the slip law and renumbered accordingly: the regularized friction
law is now eq. (12) (was 13), `V(0)` is eq. (27) (was 28), and `theta_0` is
eq. (29) (was 30). That revision did **not** change Table 1 and did **not**
resolve any of what follows.

BP8 specifies three quantities that one equation already ties together. Eq. (13)
relates `tau`, `V` and `theta`, so a case may prescribe **two** of them. The
benchmark gives all three:

| source    | quantity                                    |
|-----------|---------------------------------------------|
| Table 1   | `tau_init` = 14.6 MPa, `sigma_bar_0` = 25 MPa |
| eq. (27)  | `V(t=0)` = `V_init` = 1e-12 m/s             |
| eq. (29)  | `theta_0` = `D_RS/V_init` = 5.0e8 s         |

They are mutually inconsistent. With eq. (29)'s `theta_0`, eq. (12) gives
`f = 0.51707`, hence `tau = 12.93 MPa`, not the 14.6 MPa of Table 1. Equivalently,
holding Table 1's `tau_init` and eq. (29)'s `theta_0` forces `V(0) = 6.54e-11 m/s`,
65x the `V_init` that eq. (27) prescribes.

Three readings are possible, each sacrificing exactly one constraint:

| reading | keeps | sacrifices |
|---------|-------|------------|
| **A** literal | Table 1 `tau`, eq. (29) `theta` | eq. (27) `V` |
| **B** equilibrated | eq. (27) `V`, eq. (29) `theta`; `tau_0` -> 12.9277 MPa | Table 1 `tau_init` |
| **C** derived state | eq. (27) `V`, Table 1 `tau`; `theta_0` -> 4.02e11 s | eq. (29) `theta` |

Measured in EQquasi, dx = 50 m, `xi` = 0.2, 4000 steps:

| reading | slip at centre | peak `log10 Vmax` | `V` at first solved step |
|---------|----------------|-------------------|--------------------------|
| A (`work/c3.xi0.2`)     | 42.74 mm at 22.3 d | -6.374 @ 0.85 d | 6.542e-11 m/s |
| B (`work/bp8.equil50`)  | 23.97 mm at 23.2 d | -6.748 @ 1.64 d | 1.000e-12 m/s |
| C (`work/bp8.icfix50`)  | 37.30 mm at 22.5 d | -6.395 @ 2.37 d | 1.000e-12 m/s |

The independent whole-space BEM predicted 24.0 mm for reading B; the FEM gives
23.97 mm, agreeing to 0.1 %. Only reading B starts at `V_init` as eq. (27)
requires -- A begins 65x faster.

**Reading C is what the compsets ship**: Table 1's `tau^0` = 14.6 MPa with
`theta_0` = 4.0188e11 s derived so the fault starts at `V_init`. That is what
`taehoKim_ref` does -- its `shear_stress_2` starts at 14.6 MPa and its
`slip_rate_2` at `log10` = -12 -- and it is the closest of the three readings at
both compared stations.

The gold result is frozen in `reference/bp8/` (v1.4.7, `work/bp8.sub147`,
dx = 50 m, aging law, 30 days) and is the oracle for
`tests/e2e/test_bp8_against_gold.py`:

| quantity | ours | `taehoKim_ref` |
|----------|------|----------------|
| `tau^0` at t = 0 | 14.6000 MPa | 14.6 MPa |
| `V` at first solved step | 1.0000e-12 m/s | 1e-12 m/s |
| slip at `(0,0)`, 30 d | 37.31 mm | ~38 mm |
| slip at `(-200,0)`, 30 d | 21.99 mm | ~21 mm |
| edge/centre ratio | 0.589 | 0.553 |
| late-time pore pressure | 1.6880 MPa | ~1.70 MPa |
| peak `log10 Vmax` | -6.395 at 2.37 d | ~-6.4 |

Slip is unchanged between 23 d and 30 d, confirming that it arrests shortly
after `t_off` rather than continuing to accrue.

`shear_steady_state()` -- the BP5/BP7 helper that derives `tau^0` -- was used
briefly for reading B and has been removed from the BP8 compsets, since reading B
prescribed the one quantity Table 1 fixes. It remains in the BP5 and BP7
compsets, where it is correct.

**A dx = 50 m run is a physics test, never a submission.** Section 4.3 requires
profile nodes "with a spacing of 10 m (exactly)" from -400 to 400 m, i.e. 81
nodes; at dx = 50 m there are 17. `checkBP8Submission` rejects it with 20 errors
for that reason, which is correct -- Table 1 specifies `Delta z` = 10 m and only
the dx = 10 m run can be uploaded.

`scripts/plotICReadings.py` overlays the three (`work/ic_readings.png`). Two
things are visible there that the endpoint numbers hide:

  - **The readings differ only in the first ~4 days.** Shear stress and state at
    the centre station converge onto a common curve shortly after `t_off`; all
    the divergence is the early transient. That is consistent with the
    over-determination being the whole story for the spread between A, B and C.
  - **Our slip arrests at `t_off`; the reference's does not.** All three readings
    plateau once injection stops -- pressure diffuses away, `sigma_bar` recovers,
    the fault re-locks -- whereas `taehoKim_ref` keeps accruing slip gradually
    out to 30 d. So reading B matching ~21-24 mm is a match in *amplitude*, not
    in time history, even though the two pore pressure curves agree closely.
    That residual is **not** explained yet and should not be presented as
    agreement.

The choice changes the answer by **1.8x**, so it is not a detail. Reading A is
what the cases ship, because it takes Table 1 and section 3's eq. (29) at face
value. Reading B satisfies both equations of section 3 and adjusts only a Table 1
number, and it is the one that matches the `taehoKim_ref` comparison (~21 mm,
with the slip-rate maximum near `t_off` rather than at 0.8 d). This has been
raised with the benchmark authors; until they resolve it, `tests/
test_initial_conditions.py` pins reading A and records the size of the
inconsistency so it cannot drift unnoticed.

Note also that **Table 1 omits Poisson's ratio**, which the whole-space elastic
kernel depends on. That is a second specification gap worth raising.

### BP8-QD-S is out of scope as of 2026-08-06

Taeho Kim notified the SEAS community that the slip law has been dropped from
BP8: "In BP8, only the aging law will be used. This is to reduce the number of
simulations since we already have two variations of the fluid injection." The
updated description has been read. It removes the slip law and renumbers the
equations after it; Table 1 is unchanged, the over-determined initial condition
is unchanged, and Poisson's ratio is still absent. Both open specification
questions therefore still stand. Everything below about `friclaw = 4` is
therefore **not needed for submission**; it is kept because the two bugs it
uncovered are real, and because it bears on the open question of what the
reference actually is.

**It also undermines the conclusion recorded below.** The argument that
`taehoKim_ref` is a slip-law run rested on its state variable sitting flat at
`log10 theta` ~ 3, which the aging law cannot produce -- under the aging law
`d(theta)/dt -> 1` once slip decays and `theta` grows to ~1e6 s over 30 days, as
ours does. If the reference is an aging-law run, that flat state needs a
different explanation, and the candidates are: its state column is not `theta`
in seconds; the reference predates the scope change and was run with the slip
law; or something is wrong in that output. **This is unresolved.** Getting
Taeho's actual data files, rather than reading his curves off a comparison plot,
is the way to settle it and is the single most useful next step.

### BP8-QD-S (slip law, `friclaw = 4`) is partially working

Two defects have been fixed: the initial state is converted from `theta` to
`psi` on read (`netcdf_io.f90`), without which `exp(psi/a)` overflowed and the
first Newton solve returned NaN; and the state was being advanced over `dt`, the
fixed CFL minimum, instead of `dtev1`, the adaptive step -- about 1e5 too slow,
which froze `psi` completely.

What now works: the slip field is smooth and symmetric, and at the centre
station `psi` freezes at ~0.629 once slip decays, i.e. `theta = (Dc/V*) *
exp((psi - f*)/b)` ~ 9.1e3 s, `log10 theta` ~ 3.96. That frozen state is the
signature of the slip law, and it is close to what `taehoKim_ref` reports (flat
at `log10 theta` ~ 3) -- which was the basis for thinking the reference was QD-S.
See the scope note above: with the slip law removed from the benchmark, that
inference is no longer safe. Under the literal initial condition the slip law gives
42.78 mm at the centre against the aging law's 42.74 mm, so the evolution law
barely changes total slip -- it changes the late-time state and creep.

What does not work: `global.dat`'s `Vmax` column shows erratic excursions to
`log10 V` ~ -3.7 that are absent under the aging law. They are not Newton
failures (zero non-convergence warnings) and they leave no trace in the slip
field. `sliprate_arr` is filled at `faulting.f90:93` from the kinematic nodal
velocity difference *before* the friction solve, rather than from the solved
`v_trial` used for the moment rate, and the slip law does not damp those
transients the way the aging law does. **So for `friclaw = 4`, station and
profile output are usable and `global.dat` is not.** Not chased further because
our submission is QD-A; it must be fixed before submitting a QD-S result.

### What the discrepancy is *not*

An independent whole-space spectral boundary-integral solver, written from the
PDF and carrying no boundaries at all, reproduces the FEM under the identical
initial condition: **43.7 mm against 45.0 mm**, peak `log10 Vmax` -6.450 against
-6.489. Domain truncation therefore does **not** explain the gap, despite
section 6 of the description warning that it might. Also excluded, with
measurements: the radiation damping term (`faulting.f90:242` gives
`vs*rho/2 = 4.6244e6 Pa s/m` = `mu/(2 cs)`, applied per eq. 8), slip-magnitude
versus `slip_2` confusion, and spatial resolution (the BEM is converged to 0.5 %
across dx = 50/25/10 m).

Computational performance
---------------------

All numbers below were measured on a single shared workstation, **not** on HPC:

  - 2 x AMD EPYC 7532, 32 cores each (64 cores total), 2.4 GHz, 512 MiB L3
  - 1 TB RAM
  - Ubuntu 22.04.5, kernel 6.8.0
  - gfortran 11.4.0, Open MPI 4.1.1, netCDF 4.8.1, MUMPS 5.8.2 built locally
    via [scivision/mumps](https://github.com/scivision/mumps)
  - built with `-O3 -fopenmp`, run on **2 MPI ranks**

| compset            | `dx`   | elements | steps | solver time | wall  | peak RSS |
|--------------------|--------|----------|-------|-------------|-------|----------|
| `test.bp5.qdc`     | 4000 m | 13 500   | 102   | 15.2 s      | 16.8 s| 217 MB   |
| `test.bp7.qdc`     | 50 m   | 12 800   | 102   | 20.0 s      | 21.7 s| 233 MB   |
| `test.bp8.qdc`     | 50 m   | 12 800   | 201   | 39.5 s      | 41.2 s| 234 MB   |
| `test.bp8.qdc`     | 25 m   | 64 000   | 201   | 162.8 s     | 4m56s | 1.23 GB  |

Peak RSS is the maximum resident set size reported by `/usr/bin/time -v` for
the `mpirun` process tree, so it is roughly the per-rank peak rather than the
total across ranks.

### Parallel scaling

Measured on the 8000-element `test.bp8.qdc` case, 200 steps, `OMP_NUM_THREADS=1`,
otherwise idle machine:

| MPI ranks | wall  | s / step | speedup |
|-----------|-------|----------|---------|
| 1         | 29.4 s| 0.147    | 1.00x   |
| 2         | 25.6 s| 0.128    | 1.15x   |
| 4         | 24.1 s| 0.121    | 1.22x   |
| 8         | 23.8 s| 0.119    | 1.24x   |

Eight times the ranks buys 24 %, i.e. about 15 % parallel efficiency. At this
problem size the direct solve is dominated by work that does not distribute, so
**one rank is the sensible default**; add ranks only when the mesh is large
enough that memory, not time, is the constraint.

**Always set `OMP_NUM_THREADS`.** The build uses `-fopenmp`, so if it is unset
each MPI rank may spawn as many threads as it likes; several ranks on one node
then oversubscribe badly. `runInfo.json` records `omp_threads_per_rank = 0` when
it was unset, precisely so this is visible after the fact.

### Scaling with mesh size

Comparing two `test.bp8.qdc` runs that differ only in resolution:

| `dx`   | elements | s / step | peak RSS |
|--------|----------|----------|----------|
| 50 m   | 12 800   | 0.196    | 234 MB   |
| 25 m   | 64 000   | 0.810    | 1.23 GB  |

5x the elements costs about 4.1x the time per step and 5.2x the memory, so over
this range both are close to linear in element count. Treat any extrapolation
beyond it as an order-of-magnitude estimate only: *```EQquasi```* solves with a
direct sparse factorization (*MUMPS*), whose cost eventually grows faster than
the number of unknowns.

**Measure on an idle machine.** An earlier timing of the `dx = 25 m` case on
this same hardware gave 1549 s instead of 162.8 s -- a 10x error -- purely
because the node was under load average ~33. Check `uptime` first.

### Reading the numbers

  - "solver time" is *```EQquasi```*'s own reported time for the time loop;
    "wall" additionally includes mesh generation, matrix assembly, the initial
    *MUMPS* analysis and factorization, and netCDF output.
  - `test.bp8.qdc` runs 201 steps against the others' 102, so compare the
    per-step cost (0.15, 0.20, 0.20 s respectively), not the totals.
  - BP7 and BP8 cost slightly more per step than BP5 at a comparable element
    count because their meshes are more strongly graded along `y`.

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


## Provenance

Blessed at **v1.6.0** (`runInfo.json`). Confirmed still reproducible by the full e2e
tier on 2026-08-14, at the v1.13.0 workflow-revamp commit set.

The version gap is not staleness. These numbers are what current code
produces, checked every full-tier run; what was missing was any record of
which version produced them, so a future divergence could not be dated.
