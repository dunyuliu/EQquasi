# Compsets

A compset is one directory holding a `user_defined_params.py`, plus any static
input it ships (a rough-geometry table, a fault-geometry generator).
`create.newcase <dir> <name>` copies it into a case.

**This table is the register.** `testsys/contract/test_compsets.py` asserts it
matches the directories on disk, so it cannot drift silently. Each
`user_defined_params.py` repeats its own row in a header comment, so the status
is visible where you are editing.

## Naming (frozen 2026-08-15, user-approved)

`<benchmark>.<mode>[.<variant>].<dx_m>[.<description>]`

- **mode** — `qdc` quasi-dynamic, `fdc` fully dynamic
- **variant** — optional, open list; in use: `rough`, `kink`, `planar`, `gs`, `a`, `con`, `dip90`
- **dx_m** — on-fault element size in metres, always present
- **description** — optional free token after resolution (e.g. `leftlateral`)

Durable facts only. Gate status is deliberately *not* in the name: it changes
as the suite evolves, and a name that encodes a changing fact is a rename
waiting to happen — with `reference/` paths, `testsys/e2e/cases.py` and every
existing case to update each time.

`test.` prefixes a regression variant and is the ONLY difference from the
production grammar: where a production twin exists the test name is
`test.<production name>` exactly (`test.bp5.qdc.2000` ↔ `bp5.qdc.2000`);
where none exists the test compset still uses the full grammar. Small, cut
`nstep`, not for science.

## Register

| compset | gate | reference | published | changed |
|---|---|---|---|---|
| `bp5.qdc.2000` | e2e full (via `test.bp5.qdc.2000`) | `reference/test.bp5.qdc.2000` | Jiang et al. 2022 JGR | 2026-08-12 |
| `bp7.qdc.a.10` | e2e fast (via `test.bp7.qdc.a.10`) | `reference/test.bp7.qdc.a.10` | — | 2026-08-12 |
| `bp8.qdc.gs.10` | e2e fast (via `test.bp8.qdc.gs.10`) | `reference/test.bp8.qdc.gs.10` | CRESCENT DET (submitted) | 2026-08-12 |
| `bp1002.qdc.2500` | **e2e full, run directly** | `reference/bp1002.qdc.2500` | — | 2026-08-13 |
| `bp1002.qdc.caps.2500` | none — UNVERIFIED (live run's cycle 0 matches uncapped `bp1002.qdc.2500`) | — | — | 2026-08-15 |
| `das.qdc.10` | none | — | — | 2026-08-12 |
| `liu2020.qdc.kink.300` | none | — | Liu et al. 2020 GJI | 2026-08-14 |
| `liu2020.qdc.kink.600` | reference frozen (utilities-read; no e2e row, ~5 h/cycle) | `reference/liu2020.qdc.kink.600` | Liu et al. 2020 GJI | 2026-08-15 |
| `bp5.qdc.kink.2000` | **e2e full** | `reference/bp5.qdc.kink.2000` | — | 2026-08-15 |
| `bp1001.fdc.250` | none (ported, UNVERIFIED) | — | — | 2026-08-12 |
| `bp1001.fdc.rough.250` | none (ported, UNVERIFIED) | — | — | 2026-08-12 |
| `bp1001.qdc.rough.250` | none (ported, UNVERIFIED) | — | — | 2026-08-12 |
| `liu2020.fdc.planar.300` | none (ported, UNVERIFIED) | — | Liu et al. 2020 GJI | 2026-08-12 |
| `liu2020.fdc.rough.250` | none (ported, UNVERIFIED) | — | — | 2026-08-12 |

### Ported 2026-08-15, unverified — the five pre-`par` compset

These declared bare module-level variables (`dx = 250.0e0`) instead of the
`par` object, and `case.setup` died on `NameError: name 'par'`. Ported
mechanically on 2026-08-15 — values untouched, names prefixed, the six old
domain names mapped to `fxmin`-family — and each now passes `case.setup` and
starts the solver (planar path stepped to 40; the rough path assembled its
1.55M-node mesh from the shipped table before the smoke timeout).

**Unverified remains the operative word**: none has a `reference/`, so the
numbers are checked by nothing. The rough compset' `rough_geo_cycle.txt` was
renamed to `bFault_Rough_Geometry.txt`, the name the solver actually reads;
the old name survives only in the legacy HPC batch template.

**"gate: none" is not "broken".** `das.qdc.10` and `liu2020.qdc.kink.300` are
on the current schema and run; nothing merely checks them.

**"published"** names where results from this configuration appear. It is not
a claim that the compset reproduces that paper — `liu2020.qdc.kink.300`
demonstrably does not yet.

## Regression variants

Small and fast, for the gate, **not for science**. Outside the register.

| variant | gate | changed |
|---|---|---|
| `test.bp5.qdc.2000` | e2e fast | 2026-08-12 |
| `test.bp7.qdc.a.10` | e2e fast | 2026-08-12 |
| `test.bp8.qdc.gs.10` | e2e fast | 2026-08-12 |
| `test.bp5.qdc.dip90.2000` | e2e fast | 2026-08-15 |
| `test.stepover.qdc.1000` | e2e fast (ntotft > 1 in the fast tier) | 2026-08-15 |
| `test.stepover.qdc.con.1000` | e2e fast | 2026-08-15 |

The two-tier design of rule 7: the production compset produces the reference,
and a smaller `test.*` variant with the same `dx` but a cut `nstep` runs on
every push. `bp5.qdc.2000` and `test.bp5.qdc.2000` are not duplicates — 10000 steps
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
  `dx = 300 m`; 600 m gives 1.26 and 1200 m gives 0.63. The 64-bit-MUMPS
  blocker for 300 m is DEAD (2026-08-15): the claim was a derived estimate;
  measured, our 32-bit MUMPS factorizes 6.05e9 reals with INFOG(1)=0 in 65 s
  (MUMPS 5.x indexes the real factor array with 8-byte offsets regardless of
  intsize64). dx = 300 runs today; `work/kink300.sci` is the live attempt.
- **It is quasi-dynamic.** The paper loops EQquasi with EQdyna through EQsimu;
  every rupture in its figures 3–7 is EQdyna's. This is `par.mode = 1`,
  EQquasi alone — interseismic and nucleation only.

## Adding one

Rule 7: a production compset gets a directory and a row above; a CI variant
gets the same name with a `test.` prefix, a row in `testsys/e2e/cases.py`, and
a reference under `reference/<compset name>/` — reference directories are named
for the compset that produced them, so the three names always agree. Put the same status in the file's own header.
The contract test fails if the register and the directories disagree.
