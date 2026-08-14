# Compsets

A compset is one directory holding a `user_defined_params.py`, plus any static
input it ships (a rough-geometry table, a fault-geometry generator).
`create.newcase <dir> <name>` copies it into a case.

**This table is the register.** There is no `compsets.txt`; it was a bare list
that duplicated the directory names and said nothing useful.
`tests/contract/test_compsets.py` asserts this table matches the directories
on disk, so it cannot drift silently.

## Register

| compset | gate | reference | published | notes |
|---|---|---|---|---|
| `bp5.qdc.2000` | e2e full (via `test.bp5.qdc`) | `reference/bp5` | Jiang et al. 2022 JGR | SEAS BP5 |
| `bp7.qdc.a.10` | e2e fast (via `test.bp7.qdc`) | `reference/bp7` | — | SEAS BP7, variant a |
| `bp8.qdc.gs.10` | e2e fast (via `test.bp8.qdc`) | `reference/bp8` | submitted to CRESCENT DET | SEAS BP8, Gaussian source |
| `bp1002.qdc.2500` | **e2e full, run directly** | `reference/bp1002` | — | two-fault step-over; the only `ntotft > 1` case (rule 12) |
| `bp1001.fdc.250` | none | — | — | legacy |
| `bp1001.fdc.rough.250` | none | — | — | legacy, rough fault |
| `bp1001.qdc.rough.250` | none | — | — | legacy, rough fault |
| `das.cycle` | none | — | — | legacy |
| `liu2020.fdc.planar` | none | — | Liu et al. 2020 GJI | planar verification case of the paper |
| `liu2020.fdc.rough.250` | none | — | — | legacy, rough fault |
| `liu2020.kink.qdc` | none | — | Liu et al. 2020 GJI | 10° bend; **blocked**, see below |

**"gate: none" does not mean broken** — it means no test would notice if it
broke. `create.newcase` never validates a name; it does `os.listdir` and
`shutil.copy`.

**"published"** names where results from this configuration appear. It is not
a claim that the current compset reproduces that paper — `liu2020.kink.qdc`
demonstrably does not yet.

## Regression variants

Small and fast, for the gate, **not for science**. Deliberately outside the
register.

| variant | gate |
|---|---|
| `test.bp5.qdc`, `test.bp7.qdc`, `test.bp8.qdc` | e2e fast tier |
| `test.bp5.qdc.dip90` | none — orphaned |
| `test.stepover.qdc` | none — orphaned |

The two-tier design of rule 7: the production compset produces the reference,
and a smaller `test.*` variant with the same `dx` but a cut `nstep` runs on
every push. `bp5.qdc.2000` and `test.bp5.qdc` are not duplicates — 10000 steps
against 101.

## liu2020.kink.qdc

The EQquasi half of Liu, Duan & Luo (2020, *GJI* **220**, 598–609). A 60 × 30 km
strike-slip fault with a 10° bend, expressed through the rough-fault path
because `par.faultgeom` puts each fault on a constant-*y* plane and cannot
represent an oblique segment.

Run the geometry generator before `case.setup`:

```
create.newcase work/mykink liu2020.kink.qdc
cd work/mykink
python3 input/generateKinkGeometry.py input
./case.setup && bash run.sh
```

Set the caps to the paper's value (`par.max_norm = -100.0e6`); at the −40 MPa
default the solver refuses to start, because clamping would leave the initial
shear stress — computed at −50 MPa — no longer at steady state.

Two things block reproduction, both tracked in `PATHWAY_FORWARD.md`:

- **Resolution.** The paper requires Λ₀/dx > 2.3 and reports 2.5163 at
  `dx = 300 m`; 600 m gives 1.26 and 1200 m gives 0.63. The paper's 300 m needs
  a 64-bit-integer MUMPS (factors want ~5.5e9 reals against a 2³¹ ceiling), so
  today it runs only at resolutions its own paper rejects.
- **It is quasi-dynamic.** The paper loops EQquasi with EQdyna through EQsimu;
  every rupture in its figures 3–7 is EQdyna's. This is `par.mode = 1`,
  EQquasi alone — interseismic and nucleation only.

## Adding one

Rule 7: a production compset gets a directory and a row in the register above;
a CI variant gets `test.<benchmark>.<mode>/`, a row in `tests/e2e/cases.py`,
and a reference under `reference/<benchmark>/`. The contract test fails if the
register and the directories disagree.
