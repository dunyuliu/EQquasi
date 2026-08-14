# Compsets

A compset is one directory holding a `user_defined_params.py`, plus any static
input it ships (a rough-geometry table, a fault-geometry generator).
`create.newcase <dir> <name>` copies it into a case.

**This table is the register.** `testsys/contract/test_compsets.py` asserts it
matches the directories on disk, so it cannot drift silently. Each
`user_defined_params.py` repeats its own row in a header comment, so the status
is visible where you are editing.

## Naming

`<benchmark>.<mode>[.<variant>].<dx_m>`

- **mode** — `qdc` quasi-dynamic, `fdc` fully dynamic
- **variant** — optional: `rough`, `kink`, `planar`, `gs`, `a`
- **dx_m** — on-fault element size in metres, always present

Durable facts only. Gate status is deliberately *not* in the name: it changes
as the suite evolves, and a name that encodes a changing fact is a rename
waiting to happen — with `reference/` paths, `testsys/e2e/cases.py` and every
existing case to update each time.

`test.` prefixes a regression variant. That *is* durable: small, cut `nstep`,
not for science.

## Register

| compset | gate | reference | published | changed |
|---|---|---|---|---|
| `bp5.qdc.2000` | e2e full (via `test.bp5.qdc`) | `reference/bp5` | Jiang et al. 2022 JGR | 2026-08-12 |
| `bp7.qdc.a.10` | e2e fast (via `test.bp7.qdc`) | `reference/bp7` | — | 2026-08-12 |
| `bp8.qdc.gs.10` | e2e fast (via `test.bp8.qdc`) | `reference/bp8` | CRESCENT DET (submitted) | 2026-08-12 |
| `bp1002.qdc.2500` | **e2e full, run directly** | `reference/bp1002` | — | 2026-08-13 |
| `das.qdc.10` | none | — | — | 2026-08-12 |
| `liu2020.qdc.kink.300` | none | — | Liu et al. 2020 GJI | 2026-08-14 |
| `bp1001.fdc.250` | **TODO: broken** | — | — | 2026-08-12 |
| `bp1001.fdc.rough.250` | **TODO: broken** | — | — | 2026-08-12 |
| `bp1001.qdc.rough.250` | **TODO: broken** | — | — | 2026-08-12 |
| `liu2020.fdc.planar.300` | **TODO: broken** | — | Liu et al. 2020 GJI | 2026-08-12 |
| `liu2020.fdc.rough.250` | **TODO: broken** | — | — | 2026-08-12 |

### TODO: broken — five compsets in the pre-`par` format

`bp1001.*` and `liu2020.fdc.*` declare bare module-level variables
(`dx = 250.0e0`, `mode = 2`) instead of the `par` object every current script
reads. `case.setup` dies immediately:

```
NameError: name 'par' is not defined
```

`bp1001.fdc.250` is worse than merely old — it mixes both styles, with bare
`dx` and `par.ntotft` in the same file, so it never worked in either.

They are not ported here because there is no oracle to port them against:
none has a `reference/`, so a port would produce a plausible compset nobody
could check. Tracked in `PATHWAY_FORWARD.md`.

**"gate: none" is not "broken".** `das.qdc.10` and `liu2020.qdc.kink.300` are
on the current schema and run; nothing merely checks them.

**"published"** names where results from this configuration appear. It is not
a claim that the compset reproduces that paper — `liu2020.qdc.kink.300`
demonstrably does not yet.

## Regression variants

Small and fast, for the gate, **not for science**. Outside the register.

| variant | gate | changed |
|---|---|---|
| `test.bp5.qdc` | e2e fast | 2026-08-12 |
| `test.bp7.qdc` | e2e fast | 2026-08-12 |
| `test.bp8.qdc` | e2e fast | 2026-08-12 |
| `test.bp5.qdc.dip90` | **TODO: orphaned**, no e2e row | 2026-08-12 |
| `test.stepover.qdc` | **TODO: orphaned**, no e2e row | 2026-08-12 |

The two-tier design of rule 7: the production compset produces the reference,
and a smaller `test.*` variant with the same `dx` but a cut `nstep` runs on
every push. `bp5.qdc.2000` and `test.bp5.qdc` are not duplicates — 10000 steps
against 101.

## liu2020.qdc.kink.300

The EQquasi half of Liu, Duan & Luo (2020, *GJI* **220**, 598–609). A 60 × 30 km
strike-slip fault with a 10° bend, expressed through the rough-fault path
because `par.faultgeom` puts each fault on a constant-*y* plane and cannot
represent an oblique segment.

Run the geometry generator before `case.setup`:

```
create.newcase work/mykink liu2020.qdc.kink.300
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
  today it runs only at resolutions its own paper rejects. The directory is
  named for the paper's 300 m, not for what currently runs.
- **It is quasi-dynamic.** The paper loops EQquasi with EQdyna through EQsimu;
  every rupture in its figures 3–7 is EQdyna's. This is `par.mode = 1`,
  EQquasi alone — interseismic and nucleation only.

## Adding one

Rule 7: a production compset gets a directory and a row above; a CI variant
gets `test.<benchmark>.<mode>/`, a row in `testsys/e2e/cases.py`, and a reference
under `reference/<benchmark>/`. Put the same status in the file's own header.
The contract test fails if the register and the directories disagree.
