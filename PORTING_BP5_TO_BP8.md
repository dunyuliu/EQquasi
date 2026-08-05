# Port book: SEAS BP5 → BP8 in EQquasi

How BP8 (`BP8-QD-GS` / `-PW`, quasi-dynamic fluid injection in 3D, released
2026-07-31) was added to a code that already ran BP5 — what changed, what broke,
and what is verified versus merely working.

Read `PROJECT_RULES.md` first; this document assumes it.

---

## 1. BP8 is not a variation on BP5

They share the friction law and the solver. Almost nothing else.

| | BP5 | BP8 |
|---|---|---|
| what drives slip | far-field tectonic loading, `4e-10 m/s` | fluid injection at the fault centre |
| friction regime | velocity-weakening patch, `a=0.004 < b=0.030` | **velocity-strengthening everywhere**, `a=0.016 > b=0.010` |
| consequence | sequences of earthquakes | **entirely aseismic**, by construction |
| effective normal stress | constant, −25 MPa | **σ̄(x₂,x₃,t) = σ̄₀ − p(x₂,x₃,t)** |
| extra governing equation | none | **2-D pore pressure diffusion on the fault** |
| outside the frictional zone | creeps at `1e-9 m/s` | **no slip**, zero fluid flux |
| fault half-length | 60 km | **400 m** |
| cell size | 2000 m | **10 m** |
| duration | earthquake cycles, years | **30 days**, one run |
| termination | slip-rate threshold | **simulation time** |

The one genuinely new piece of physics:

```
∂p/∂t = α ∇²p + q_inj(t)/(βφ) · δ(x₂)δ(x₃),      α = k/(φβη) = 0.05 m²/s
```

solved on the fault plane, fed back as `σ̄ = σ̄₀ − p`.

## 2. The leverage: start from BP7, not BP5

`bp == 7` already existed, and its **geometry is BP8's**: ±400 m frictional
region in a ±500 m box, 10 m cells, 25 MPa, ρ = 2670, cₛ = 3464, Dc = 0.5 mm.
`library_output.f90` already switched station filenames to metres for `bp == 7`.

BP7 lacks only the fluid: it nucleates with an artificial shear-stress bump
(`nucleation1`, `faulting.f90:502`) on a velocity-*weakening* disc.

So the port is: **reuse BP7's mesh, stations and output plumbing; replace
artificial nucleation with pore-pressure forcing; flip a/b to velocity
strengthening.**

**The catch, and the main lesson of this port: BP7 could not run.** See §7.
Inheriting a code path also inherits whatever has never executed in it.

## 3. Three contracts you must not break

Each now has a guard in `tests/test_io_contracts.py`.

### 3.1 `model.txt` is positional

`case.setup:create_model_input_file()` writes N unlabelled lines;
`read_input.f90:readmodel` reads them in exactly that order. No keys. Insert a
line in one place only and every later parameter shifts — the run does not
crash, it produces wrong physics.

**Append only, never insert** (rule 4). BP8 appends six lines and reads them
with `iostat`, so an older `model.txt` still loads and leaves `fluid_src = 0`:

```fortran
read(1002,*,iostat=ios) fluid_src
if (ios == 0) then
    read(1002,*) fluid_q0, fluid_toff, fluid_tend
    ...
endif
```

### 3.2 On-fault variables travel by name, land by magic index

`case.setup` writes named netCDF variables; `netcdf_read_on_fault` reads them by
name into `fric(<index>, node, fault)`:

| idx | quantity | idx | quantity |
|---|---|---|---|
| 7 | initial normal stress | 13 | `r0` |
| 8 | initial shear stress | 20 | state variable |
| 9 | `a` | 46 | initial slip rate |
| 10 | `b` | **6** | **pore pressure** (was declared, unused, hardwired to 0) |
| 11 | `Dc` | **51, 52** | **Darcy velocity, strike / dip** (new) |
| 12 | `v0` | | |

Node ordering is `n = (ix-1)*nzt + iz`. The pore-pressure solver relies on this
to build its 2-D stencil without a separate grid.

That `fric(6,...)` was *already* wired into the creeping branch as
`tnrm0 = tnrm + fric(6,i,ift)` is why effective-stress coupling is one line:

```fortran
if (bp == 8) tnrm0 = min(tnrm0 + fric(6,i,ift), 0.0d0)
```

### 3.3 Fixed-width output arrays indexed by magic number

`fltsta(<n>, nstep, nonmx)` and `globaldat(10, nstep)` are indexed by bare
integers. BP8 needed four more station fields, so `fltsta` went 10 → 14.
Nothing warns you if you index past the allocation. Now a test does.

## 4. What changed

All additive and guarded on `bp == 8` (rule 6); BP5/BP7 paths byte-identical.

| file | change |
|---|---|
| `src/porepressure.f90` | **new** — `pore_pressure_init` / `pore_pressure_update` |
| `src/globalvar.f90` | fluid parameters, `pf`, `pfActive`, `pfWt`, `dtmax` |
| `src/read_input.f90` | six appended `model.txt` records, read with `iostat` |
| `src/faulting.f90` | effective-stress coupling; four extra `fltsta` fields |
| `src/solveTimeLoopMUMPS.f90` | init + update calls; `dtmax` cap; time-based exit |
| `src/library_output.f90` | 11-field station files, BP8 `global.dat`, signed-metre names |
| `src/eqquasi.f90` | `fltsta` 10 → 14 |
| `scripts/defaultParameters.py` | fluid parameters, all defaulting to off |
| `scripts/case.setup` | append the six lines |
| `case_input/bp8.qdc.gs.10/` | production compset, 10 m |
| `case_input/test.bp8.qdc/` | smoke compset, 50 m |

Plus one **shared-code bug fix** that was not BP8-specific: `src/meshgen.f90`.
See §7.

## 5. Numerics decisions

**Time step.** BP5's `dtev1 = ξ·Dc/V_max` runs away in BP8: ~10⁷ s at
`V_init = 1e-12`. BP8 caps it at `dtmax = 500 s`. Not arbitrary — that is also
the explicit-diffusion stability limit `Δx²/(4α)` at `Δx = 10 m`, and the max
time step in the benchmark's own worked example.

**Diffusion scheme.** Explicit FTCS with sub-stepping to hold
`dt_sub ≤ 0.2·Δx²/α`. Chosen for auditability; the diffusion solve is
negligible beside the MUMPS factorisation at these sizes.

**A discrepancy in the benchmark description.** Eq. (20) writes the Gaussian
source as `(1/L²)·exp(−r²/2L²)`, which integrates to `2π`, not 1 — inconsistent
with the closed form in eq. (22), which is the convolution of the point solution
with a *unit-normalised* Gaussian. The implementation follows eq. (22)/(25) and
normalises the weights **discretely** so `Σw = 1` exactly, removing both the
analytic ambiguity and the discretisation error of the smeared source.

**`r_well` is not in Table 1.** "Typically 1–4 inches", never fixed, though the
Peaceman well index depends on it logarithmically. Defaulted to 0.0508 m and
exposed as `par.fluid_rwell`.

**No-slip outside Ω_f.** The creeping branch divides by `load_slip_rate`, so a
literal zero gives `log(0)`. BP8 uses the benchmark's own `V_zero = 1e-20`.

## 6. Verification

Three independent angles, because any one alone is weak.

**1. Convergence to the closed form (eq. 22).** Error in pore pressure at the
injection point:

| t (s) | `dx`=50 m | `dx`=25 m |
|---|---|---|
| 25 000 | 5.11 % | 1.51 % |
| 100 000 | 5.60 % | 1.37 % |

Halving `dx` cuts error ~3.9× — **second order**. The residual is source
under-resolution: `L_gauss = 50 m` is one cell at `dx = 50 m`.

**2. Discrete mass conservation — independent of the Green's function.** After
shutoff the zero-flux boundary traps all fluid in Ω_f, so `p` must relax to a
uniform value fixed only by injected volume. Predicted 1.4948 MPa, measured
**1.5013 MPa** at 30 days: **+0.43 %**. This exercises the operator, the source
normalisation and the boundary condition together.

**3. Friction solve against hand derivation.** With `τ⁰ = 14.6 MPa`,
`σ̄ = 25 MPa`, `θ = Dc/V_init`, the regularised law gives `V = 6.54e-11 m/s`;
the code reports **6.55e-11**.

> **This check was hollow — the trap to avoid when porting a case file.** It
> verified the solver against an input that was itself wrong. `τ⁰`, `V_init`
> and `θ⁰` are not independent; eq. (9) ties them together, so a case may
> prescribe **two**. The BP5 case sets `θ⁰ = Dc/V_init` and *computes* `τ⁰`
> from it (`shear_steady_state()`, `case_input/bp5.qdc.2000/
> user_defined_params.py:103-108`). Porting to BP8 we kept BP5's `θ⁰` line and
> swapped its computed `τ⁰` for Table 1's 14.6 MPa. Each line looks right in
> isolation; together they over-determine the state and the solver silently
> started at `V = 6.5e-11 m/s`, 65× `V_init`, on a fault 1.67 MPa weaker than
> specified. Reproducing a number you derived from the same wrong inputs proves
> only that the arithmetic is consistent. Fixed in v1.4.6; guarded by
> `tests/test_initial_conditions.py`.

**Known bias, first order in `dx`.** Check 2 matches the *discrete* area
`(17·50)² = 722 500 m²`, not the continuum `800² = 640 000 m²`. Pressure lives
at nodes and every node carries a full `dx²` cell including boundary nodes, so
the storage area overshoots by one cell per side: 13 % area, hence **11 % low**
on the asymptotic pressure at `dx = 50 m`, ~2.5 % at `dx = 10 m`. Not covered by
the second-order result above, which is measured before the front reaches the
boundary.

**Not verified.** Slip rate peaks ~4.5 days *after* injection stops. That is the
expected signature of a slip front that keeps propagating while the growing
slipping patch loads the centre elastically — but BP8 was released in July 2026
and no community comparison exists. Expected is not verified.

## 7. What the port exposed in existing code

The reusable part of this document.

1. **`meshgen.f90:500` — far-field tagging clobbered fault elements.** Elements
   near the fault are tagged `et=2`, then elements within `2·dymax` of
   `ymin/ymax` are tagged `et=3`. The second test ran unconditionally. Only
   `et=2` elements give fault nodes their mass, so when the far-field band
   covered the fault band, **every** fault node went massless and the run
   aborted naming *one corner node* — making a total failure look local.
   `dymax = min(12·dx, 3 km)`, so a ±500 m model at `dx = 50 m` gets a 1200 m
   band against a 500 m half-width. BP5 never sees it (`dymax` saturates at
   3 km against ±50 km). **This blocked both BP7 and BP8.**
2. **BP7 could not run.** `case_input/bp7.qdc.a.10` set `par.xmin/...` but
   `case.setup` reads `par.fxmin/...`. The assignment silently created unused
   attributes, so BP7 inherited BP5's 120×100×60 km domain — at `dx = 10 m`.
   Also affected `bp5.qdc.2000` (masked, values coincide with defaults) and
   `das.cycle`.
3. **BP7 had no smoke case, no oracle, no CI entry** — which is why (2)
   survived. The only BP7 compset was a 35-hour, 20-CPU production run.
4. **The release gate could not fail.** `check.test.py` printed `FAIL ` while CI
   grepped for `error|failed|fatal|abort` — no match — all under `|| true`.
5. **Five compsets** (`bp1001.*`, `liu2020.*`) do not use the `par` object at
   all and fail on import. Predates the current convention. Still open.
6. `install.eqquasi.sh`'s `local` branch was dead code (two malformed bash
   tests); the `local` makefile target linked a system ScaLAPACK that cannot
   exist on a host needing a local MUMPS build.

**The generalisable rule: before porting *from* a benchmark path, verify that
path actually runs.** Inheriting BP7's plumbing meant inheriting a code path
that had never executed. Two of the six items above were only discovered because
BP8 forced BP7 to run for the first time.

## 8. Process notes

Worth recording because they cost real time.

- **Fix the gate before trusting it.** Items 3 and 4 meant the project's safety
  net was decorative. That work came first, and it immediately caught two
  defects in the BP8 change itself (an unregistered compset, and a test whose
  verdict depended on leftover state from a previous run).
- **Never pipe a long-running measurement through a filter you truncate.**
  `mpirun ... | grep x | head -1` sends SIGPIPE and kills the run; `/usr/bin/time`
  then cheerfully reports the resource usage of a truncated process. It produced
  a plausible-looking 1.6 s / 46 MB that was pure fiction.
- **Measure on an idle machine.** The `dx = 25 m` case timed at 1549 s under
  load average ~33 and 162.8 s idle — a 10× error that produced a wrong scaling
  law in the README before it was caught and retracted.
- **One run, one directory.** A relaunch script that `rm -rf`'d its working
  directory destroyed a still-running 30-day job that had 65 minutes left on its
  timeout.
