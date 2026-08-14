# EQquasi Project Rules

This book has two tiers. **Part I** states the general principles — the ones
that would hold for any numerical-modeling project with frozen oracles,
shared compute, and parallel agents, stated so they'd make sense to someone
who has never seen this repo. **Part II** is where each principle actually
binds in EQquasi: the specific file, the specific incident, the specific
command that checks it. Read Part I once; open a Part II rule only when it's
load-bearing.

## Index

**Part I — General principles**

- [G1. Do the smallest amount of work that solves the problem](#g1-do-the-smallest-amount-of-work-that-solves-the-problem)
- [G2. A guard stops and explains; it never degrades silently](#g2-a-guard-stops-and-explains-it-never-degrades-silently)
- [G3. A gate is only real if it runs, fails loud, and is evaluated at the right configuration](#g3-a-gate-is-only-real-if-it-runs-fails-loud-and-is-evaluated-at-the-right-configuration)
- [G4. A positional contract between a writer and a reader can only be extended by appending](#g4-a-positional-contract-between-a-writer-and-a-reader-can-only-be-extended-by-appending)
- [G5. A value duplicated across independent sources of truth needs one check before every new entry](#g5-a-value-duplicated-across-independent-sources-of-truth-needs-one-check-before-every-new-entry)
- [G6. New behavior is an additive branch; the shared/default path is never edited to fit it](#g6-new-behavior-is-an-additive-branch-the-shareddefault-path-is-never-edited-to-fit-it)
- [G7. A configuration doesn't exist to the safety net until it's registered where the net looks](#g7-a-configuration-doesnt-exist-to-the-safety-net-until-its-registered-where-the-net-looks)
- [G8. Reference data is ground truth: read-only, only grows, and dead if nothing reads it](#g8-reference-data-is-ground-truth-read-only-only-grows-and-dead-if-nothing-reads-it)
- [G9. A build/environment claim is a hypothesis until it's run on the target host](#g9-a-buildenvironment-claim-is-a-hypothesis-until-its-run-on-the-target-host)
- [G10. Docs move with the code, in the same change](#g10-docs-move-with-the-code-in-the-same-change)
- [G11. A comparison is only as trustworthy as its controls](#g11-a-comparison-is-only-as-trustworthy-as-its-controls)
- [G12. A codebase-specific footgun gets documented once, with a mechanical check where possible](#g12-a-codebase-specific-footgun-gets-documented-once-with-a-mechanical-check-where-possible)
- [G13. In a shared checkout, isolate your work](#g13-in-a-shared-checkout-isolate-your-work)
- [G14. A shared compute resource is not yours alone](#g14-a-shared-compute-resource-is-not-yours-alone)
- [G15. A finding is a hypothesis until checked against the source](#g15-a-finding-is-a-hypothesis-until-checked-against-the-source)

**Part II — Project-specific rules**

1. [Minimal changes; no new files until necessary](#1-minimal-changes-no-new-files-until-necessary) (G1)
2. [No silent fallbacks, swallowed errors, or placeholder data](#2-no-silent-fallbacks-swallowed-errors-or-placeholder-data) (G2)
3. [The regression gate is the release gate — four benchmarks, real exit codes, the right configuration](#3-the-regression-gate-is-the-release-gate--four-benchmarks-real-exit-codes-the-right-configuration) (G3)
4. [model.txt is a positional contract — never insert, always append](#4-modeltxt-is-a-positional-contract--never-insert-always-append) (G4)
5. [fric()/on_fault_vars magic indices have one authoritative source per direction](#5-fricon_fault_vars-magic-indices-have-one-authoritative-source-per-direction) (G5)
6. [New bp values are additive branches, never edits to existing bp branches](#6-new-bp-values-are-additive-branches-never-edits-to-existing-bp-branches) (G6)
7. [New compsets: naming, placement, and registration](#7-new-compsets-naming-placement-and-registration) (G7)
12. [The regression gate must include a multi-fault case](#12-the-regression-gate-must-include-a-multi-fault-case) (G7)
8. [Reference results are read-only, only grow, and every file has a reader](#8-reference-results-are-read-only-only-grow-and-every-file-has-a-reader) (G8)
9. [Build/environment claims must be verified on the target host, not assumed](#9-buildenvironment-claims-must-be-verified-on-the-target-host-not-assumed) (G9)
10. [Docs move with the code, in the same change](#10-docs-move-with-the-code-in-the-same-change) (G10)
11. [Compare like with like: pin the version, vary one thing, stress the assumption](#11-compare-like-with-like-pin-the-version-vary-one-thing-stress-the-assumption) (G11)
13. [Known Fortran landmines in this codebase must not recur](#13-known-fortran-landmines-in-this-codebase-must-not-recur) (G12, and 13a also instantiates G6)
14. [Shared-checkout discipline: worktrees, not branches; explicit paths, not `-A`](#14-shared-checkout-discipline-worktrees-not-branches-explicit-paths-not--a) (G13)
15. [Operational safety on the shared compute host](#15-operational-safety-on-the-shared-compute-host) (G14)
16. [A subagent's or audit's finding is a hypothesis until checked against source](#16-a-subagents-or-audits-finding-is-a-hypothesis-until-checked-against-source) (G15)

Rule numbers are cited in commits, reports, and the test suite itself
(`tests/contract/test_reference_gold_is_referenced.py` cites 8a,
`tests/unit/test_physical_invariants.py` cites 12,
`tests/regression/test_code_convention_landmines.py` cites 13a/13b/13c) —
they never change even when a rule's position in this document does. This
file was reorganized general-first/specific-second on 2026-08-12; every rule
number and every incident is unchanged from before that pass, only the
ordering and framing changed. See Rule-book health at the bottom for what
else changed that day.

---

# Part I — General principles

These sixteen rules exist because eleven distinct project-specific incidents
happened here, but every one of them is a specific case of a principle that
would hold on a different codebase with the same shape: frozen numerical
oracles, more than one contributor, a shared machine. State the principle
first, because the specific incident is what taught it to us, not what it is.

### G1. Do the smallest amount of work that solves the problem

Fold new logic into the file that already owns the concern; don't create a
new file, and don't refactor unrelated code, in the same change that does
something else. Every other principle in this book bounds what work you do;
this one bounds how much of it you do at all — the same change that adds a
feature and reorganizes a shared file becomes impossible to bisect against
whatever regression gate is watching (G3).

→ Rule 1.

### G2. A guard stops and explains; it never degrades silently

Missing input, a missing binary, a missing or malformed config: all of these
fail loudly, at the point they're discovered. No substituted default, no
skipped step, no partial result that looks like a complete one. Silent
degradation is not a lesser violation of this than a crash — it's worse,
because a crash is at least visible.

→ Rule 2.

### G3. A gate is only real if it runs, fails loud, and is evaluated at the right configuration

A merge or release gate is only as real as three things holding
simultaneously: it actually executes end to end (not skipped, not swallowed
by a build failure upstream); its pass/fail signal is read from its real exit
code — never a piped command's exit code, never a keyword grep that can miss
the tool's actual output string; and it is run at the configuration it
claims to certify, not a faster stand-in with different parameters that
happens to share a name.

→ Rule 3.

### G4. A positional contract between a writer and a reader can only be extended by appending

When one file writes a data format and another reads it back, and the two
share no schema — fields are matched by position, not by name — that
contract can only be safely extended by appending a new field at the end.
Inserting a field anywhere else silently shifts every field after it, and
because there are no names to check against, both sides "succeed": no crash,
just wrong values in the wrong variables.

→ Rule 4.

### G5. A value duplicated across independent sources of truth needs one check before every new entry

Some values have no single schema to own them and end up written in more
than one independent place by convention rather than by generation. Keeping
those copies in sync means: grep for what's already used before allocating a
new one, write every copy in the same change, and never repurpose an
existing entry for a new meaning — that would silently corrupt every
consumer still relying on the old one.

→ Rule 5.

### G6. New behavior is an additive branch; the shared/default path is never edited to fit it

New behavior is added as a new branch alongside the existing ones. The
default (or degenerate-case) path that every existing case already exercises
is never modified to accommodate something new, and a special case is never
carved out for what used to be the only case. One code path, with the
simple/default configuration as its degenerate instance, is what makes a
passing regression suite *proof* that old and new behavior are independent —
not proof by coincidence.

→ Rules 6 and 13a.

### G7. A configuration doesn't exist to the safety net until it's registered where the net looks

Adding a new case to the project (a new benchmark, a new compset, a new
variant) does nothing for correctness by itself — it has to be registered in
every place the safety net actually looks: the production registry and the
fast CI-regression registry both. And the registry itself has an obligation
back: the CI-regression tier must include at least one case that actually
stresses the project's structurally hard assumptions, not only the cases
that are cheapest to run. A gate that only ever exercises the easy case is a
gate with a permanent blind spot.

→ Rules 7 and 12.

### G8. Reference data is ground truth: read-only, only grows, and dead if nothing reads it

Frozen oracles are the thing every other check is measured against, so
nothing writes through them — not a "just this once" regeneration, not a
side effect of unrelated work. And the corollary that's easy to forget: a
reference artifact that no test actually reads is not an asset sitting in
reserve, it's dead weight that *looks* like a safety net until someone
checks whether anything consults it.

→ Rule 8, including 8a.

### G9. A build/environment claim is a hypothesis until it's run on the target host

"This should work on host X" is a hypothesis, not a fact, until it has
actually been run there once. A build should fail immediately and legibly on
an unmet precondition — an unset variable, a missing library, a wrapper
binary that doesn't exist on this host's flavor of Linux — rather than
silently guessing a default and proceeding.

→ Rule 9.

### G10. Docs move with the code, in the same change

Documentation is not correct until it matches the code, and "matches the
code" has an expiration date the instant the code changes again. Docs move
in the same commit as the code that makes them true — never queued as a
follow-up, because the follow-up is where every doc-vs-code drift in this
project's history was born.

→ Rule 10.

### G11. A comparison is only as trustworthy as its controls

A comparison between two runs, two configurations, or two versions is only
as trustworthy as three controls held simultaneously: the two results come
from the same binary/version (not an old one and a new one that happen to
look similar), exactly one variable differs between them, and the
configuration under test is chosen because it would actually expose a
violated assumption — not because it happens to satisfy that assumption by
coincidence and therefore hides the bug that only shows up when it doesn't.

→ Rule 11.

### G12. A codebase-specific footgun gets documented once, with a mechanical check where possible

Every nontrivial codebase accumulates a handful of landmines that are
specific to it — a shadowed name, a reserved unit number, a flag with a
surprising default — each one capable of costing a full debugging session to
whoever finds it "the hard way." The fix is not vigilance, it's a single
place these are written down, with a grep-able mechanical check wired in
wherever the landmine is expressible as a pattern.

→ Rule 13.

### G13. In a shared checkout, isolate your work

When more than one contributor or agent shares one checkout, isolation is
the first obligation: your own line of work lives in your own worktree, not
a branch switch in the shared tree, and you stage exactly the paths your
change touches, never a blanket operation that could sweep up someone else's
unfinished, unreviewed work into your commit.

→ Rule 14.

### G14. A shared compute resource is not yours alone

On a machine other people and other agents are actively using, default to
the smallest footprint that does the job — the fewest ranks, the fewest
threads, artifacts confined to one place — and know the specific operational
failure modes that look correct but aren't, because those are exactly the
ones that take down someone else's run without warning.

→ Rule 15.

### G15. A finding is a hypothesis until checked against the source

A finding — from a subagent, an automated audit, a teammate's report, or a
rule book like this one — is evidence to go check, not a conclusion to act
on. It gets verified against the primary source (the paper, the equation,
the file and line) before it drives a fix, a routed report, or a change to
this file.

→ Rule 16.

---

# Part II — Project-specific rules

Each rule below states, in one line, which Part I principle it instantiates,
then gets specific: the file, the incident, the command that checks it.

## 1. Minimal changes; no new files until necessary

**Instantiates**: G1.

Smallest edit that solves the problem. Fold new logic into the file it
belongs with (e.g., a new `bp==8` block goes inside the existing `if (bp ==
N)` chain in `src/faulting.f90`, `src/library_output.f90`,
`src/solveTimeLoopMUMPS.f90` — not a new source file). Never refactor
unrelated code in the same change.

**Incident**: a pore-pressure-diffusion feature that also reorganized
`globalvar.f90` made the diff impossible to bisect against the regression
gate (rule 3).

The same applies to names, not just files. A new variable, parameter or
constant is a new thing to maintain and to keep consistent, so it gets the
same scrutiny: reuse an existing one, or write the value inline, before
inventing a name. In particular `par` carries a defined schema
(`scripts/defaultParameters.py`); a key that is not in that schema is not a
parameter, however well it reads. Compset-local geometry belongs at module
level, as `FAULT_A`/`FAULT_B` do in `case_input/bp1002.qdc.2500`.

**Incident**: the Liu et al. (2020) kink compset was written with
`par.kinkAngleDeg` and `par.kinkX` hung off `par` and several new
module-level constants, none of them flagged. They worked only because the
generator read them with `getattr(..., default)` -- the schema never knew
about them.

**How to apply**: before adding a file, grep for whether an existing file
already owns that concern (`src/faulting.f90` owns fault-node physics,
`src/library_output.f90` owns benchmark-format output, `scripts/lib.py` owns
shared Python helpers). Before adding a file OR a name, say so and get
agreement first; do not introduce either as a side effect of doing something
else.

**Tier**: judgment (2).

---

## 2. No silent fallbacks, swallowed errors, or placeholder data

**Instantiates**: G2.

Missing input, binary, or config fails loudly. No substituted default, no
skipped step, no warn-and-proceed.

**Fixed since seeding (verified 2026-08-12)**: `install.eqquasi.sh`'s
malformed `[ $MACHINE == "local"]` / `[! -f ...]` bash tests and the
`.github/workflows/test.yml` `python3 `pytest -m e2e` >> testRunLog.txt || true`
that discarded `pytest -m e2e`'s exit code are both gone from the current files;
`tests/contract/test_release_gate.py` now guards the exit-code/grep path
directly (see rule 3). Left here as the incident record — the failure mode (a
syntax slip or a swallowed exit code makes a broken branch look green) is
what future edits to either file must not reintroduce.

**2a. Silent degradation is the same violation, not a lesser one.** A
multi-fault mesh where a fault plane fell between mesh lines got zero nodes
silently — verified at dx=1000 with a 5 km fault offset (5 cells exactly, an
assumption-satisfying configuration, see rule 11), it worked; at dx=2000 the
offset fell between mesh lines and one fault got zero nodes. The run
completed and wrote output as if nothing were wrong, attributing 11 m of slip
to what was effectively a phantom single fault. Every place that counts
fault nodes per fault (`nonfs(:)` in `src/read_input.f90`, the per-fault loops
in `src/mesh4num.f90`/`src/meshgen.f90`) must `stop` with the fault index and
node count if any fault ends up with zero nodes — never proceed.
**Not currently guarded**: `grep -n "stop" src/read_input.f90
src/mesh4num.f90` shows no check on `nonfs(i) == 0` today (verified
2026-08-12) — this is an open gap, not a hypothetical; see Violations.

**How to apply**: a new pore-pressure input file, a new `bp==8` compset
parameter, a new solver path, or a new per-fault array must `stop` with a
descriptive message on missing/malformed/empty input — following the
existing pattern at `src/read_input.f90:23-34` (`INQUIRE` + explicit `stop`)
— never fall back to a default value, an empty array, or a partial result.

**Tier**: mechanical for the two now-fixed syntax bugs (`bash -n
install.eqquasi.sh`, `grep -n '|| true' .github/workflows/test.yml`) and,
once it exists, for the zero-node guard (`grep -n "nonfs(i) == 0" src/*.f90`);
judgment for new code in general.

---

## 3. The regression gate is the release gate — four benchmarks, real exit codes, the right configuration

**Instantiates**: G3.

`reference/bp5/`, `reference/bp5.dip90/`, `reference/bp7/` and
`reference/bp8/` are the project's safety net (moved here from
`test.reference.results/` when the reference layout was unified, git
`ab05778`/`a771065` — that older path no longer exists; do not look for it).
Before merging any change to `src/*.f90` or to
`scripts/defaultParameters.py`/`scripts/case.setup`:

```
python3 -m pytest tests/            # unit + contract + regression, ~4 s
python3 -m pytest -m e2e tests/     # builds and runs BP5/BP5-dip90/BP7/BP8
python3 scripts/plotAgainstGold.py bp5      <rundir>
python3 scripts/plotAgainstGold.py bp5.dip90 <rundir>
python3 scripts/plotAgainstGold.py bp7      <rundir>
python3 scripts/plotAgainstGold.py bp8      <rundir>
```

All four `plotAgainstGold.py` invocations must print `matches gold` and exit
0, and every `pytest` invocation must exit 0. Read the real exit status
(`$status` in tcsh, `$?` in bash) — never a piped command's. `python3 -m
pytest -q | tail` reports `tail`'s exit code, not pytest's; a red suite piped
through `tail` looks identical to a green one. This is the same failure mode
CI already hit once (below); do not reintroduce it in a local workflow.

**Fixed since seeding (verified 2026-08-12)**: the two gaps originally
recorded here — `tests/e2e/test_benchmarks.py`'s literal `'FAIL '` not matching CI's grep
pattern, and `pytest -m e2e`'s exit code being discarded by `|| true` — are both
gone from `.github/workflows/test.yml` (now `set -o pipefail`, an explicit
`grep -qE 'FAIL'` check, and a `CHECK SUMMARY` completion check), and are
themselves guarded by `tests/contract/test_release_gate.py`. `tests/e2e/test_benchmarks.py`
owns every benchmark's comparison through `tests/e2e/cases.py`'s
`compare_series`/`compare_netcdf`, at `rtol = 1e-9` with a floor scaled to
each quantity's own vector (`floor_rtol = 1e-9`). The earlier text here
described `compare_nc_files`/`compare_txt_files` at `rtol=atol=1e-3`; those
functions no longer exist and the figure was six orders of magnitude off the
gate as it now runs, so anyone sizing a change against the rule book was
reasoning from a number the code had abandoned.
`reference/bp8/summary.json` plus `tests/e2e/test_bp8_against_gold.py`
own BP8's, at tolerances near round-off (`RTOL_SLIP = 1e-4`, exact for the
fault-plane snapshot).

**BP8 must be gated at its gold configuration**, not the fast regression
compset's defaults: `xi = 0.2`, run to `fluid_tend` rather than an `nstep`
cap — the gold run took 5301 steps ≈ 30.0 days
(`reference/bp8/runInfo.json`: `"steps_completed": 5301`).
`case_input/test.bp8.qdc/user_defined_params.py` defaults to `nstep = 200`
(~1.14 days) for the *fast* regression tier and must never be compared
directly to the BP8 gold — this has already misled two agents.
`tests/e2e/test_bp8_against_gold.py` gets this right today (`run_dir`
fixture patches `xi`/`nstep`/`nt_out` before running `test.bp8.qdc`); any new
script that runs BP8 against the gold must copy that pattern, not
`test.bp8.qdc`'s checked-in defaults as-is.

**Route**: a regression in one of the four benchmarks → `lars-eriksson`
(code) or `haruto-nakamura` (if it's the gate itself that's wrong).

**Tier**: mechanical — run the commands above and read stdout/exit code.

---

## 4. model.txt is a positional contract — never insert, always append

**Instantiates**: G4.

`scripts/case.setup`'s `create_model_input_file` writes fixed-order lines;
`src/read_input.f90`'s `readmodel` reads them back with unlabeled
`read(1002,*)` statements in the same order, matched positionally, not by
name. `stations.txt` (`create_station_input_file` / `readstations1`/`2`) and
`fric.txt` (`readfric`, format branches on `friclaw`) follow the same
convention.

Any new parameter (e.g., a pore-pressure diffusivity, an injection rate, a
diffusion-grid spacing) **must be appended after the last existing
`f.write(...)`/`read(1002,*)` pair on both sides**, never inserted in the
middle — an insertion silently shifts every field after it and both sides
will "succeed" (no crash) while reading garbage into the wrong variable.

**Incident risk**: this is exactly the kind of change that passes a naive
smoke test (`case.setup` runs, `eqquasi` starts) while producing wrong
physics, because Fortran's list-directed `read(unit,*)` has no field names to
catch a shift.

**How to apply**: when adding a line, edit both files in the same change, add
it as the new last line, and rerun the regression gate (rule 3) to confirm
the *unrelated* existing fields still parse correctly.

**Tier**: judgment (2) — requires reading the diff of both files side by
side; not mechanically checkable without a schema, which does not exist (see
Unenforceable rules).

---

## 5. fric()/on_fault_vars magic indices have one authoritative source per direction

**Instantiates**: G5.

`fric(1:100, node, fault)` is indexed by bare integers with no named
constants. The same indices are duplicated in three places that must never
disagree: `scripts/defaultParameters.py` (Python-side array),
`scripts/case.setup` (`netcdf_write_on_fault_vars`, both the `createVariable`
call and the `var[:,:] = par.on_fault_vars[:,:,N]` assignment), and
`src/netcdf_io.f90` (the corresponding Fortran read).

- **Pick an index not already in use.** Grep `fric(<N>,` across `src/*.f90`
  first.
- **Write it in exactly three places in the same change.**
- **Never repurpose an existing index** for a new meaning — it will silently
  corrupt every compset that still writes to it under the old meaning.

**Route**: an addition that misses one of the three files is a doc-vs-code
drift bug → `sophia-okafor`; a mismatched index value across the three (data
corruption, not a doc issue) → `lars-eriksson`.

**Tier**: judgment (2) for new additions; mechanical only in the narrow sense
that `grep -n "fric(<N>," src/*.f90` before assigning a new index is a
checkable precondition.

---

## 6. New bp values are additive branches, never edits to existing bp branches

**Instantiates**: G6.

Existing benchmark selection follows one pattern throughout `src/*.f90`:

```fortran
if (bp == 5) then
    ...
endif
if (bp == 7) then
    ...
endif
```

A new benchmark adds a new `if (bp == N) then ... endif` block adjacent to
the existing ones in each file that branches on `bp` — **never modify the
body of an existing `bp == 5`/`bp == 7`/`bp == 8` block** to accommodate new
logic, and never change the default (no-`bp`-match) code path.

**Rationale**: this is what makes rule 3's gate meaningful — if new-benchmark
work can only add branches, a passing regression run after that change is
proof the two are actually independent, not proof by coincidence. Rule 13a
makes the identical argument one level down, for `ntotft` instead of `bp`.

**How to apply**: `git diff` before committing should show 0 modified lines
inside any pre-existing `if (bp == N)` block in `src/faulting.f90`,
`src/library_output.f90`, `src/solveTimeLoopMUMPS.f90`.

**Tier**: mechanical in principle (a script rejecting any `-`/`+` line pair
inside an existing `if (bp == N)` block's byte range); not implemented today
— see Unenforceable rules.

---

## 7. New compsets: naming, placement, and registration

**Instantiates**: G7, together with rule 12 below.

- **Production compset** → new directory `case_input/<benchmark>.<mode>.<res>/`
  (matching `bp5.qdc.2000`, `bp7.qdc.a.10`, `bp8.qdc.gs.10`), with an entry
  appended to `case_input/compsets.txt`.
- **CI regression compset** → a second, smaller/faster directory named
  `test.<benchmark>.<mode>` (mirroring `test.bp5.qdc` vs `bp5.qdc.2000`: cut
  `nstep` and coarsen `dx`), added to `nameList`/`coreNumList` in
  `tests/e2e/cases.py` (or wired into a `tests/e2e/*.py` regression, for
  benchmarks that moved off `tests/e2e/cases.py` — see rule 12), with its
  reference output at `reference/<benchmark>/`.

`scripts/create.newcase` does **not** validate `compset` against
`compsets.txt` — it just does `os.listdir` + `shutil.copy` from
`case_input/<compset>/`. `compsets.txt` is documentation-only, not an
enforced allowlist; keep it accurate rather than relying on it to reject
typos.

**Route**: `compsets.txt` drifting from the actual `case_input/` directories
present → `sophia-okafor`.

**Tier**: mechanical (`diff <(ls case_input | grep -vE
'compsets.txt|user_defined_params.py|^test\.') <(sort case_input/compsets.txt)`)
— not implemented today, see Unenforceable rules.

---

## 12. The regression gate must include a multi-fault case

**Instantiates**: G7, together with rule 7 above — a registry that never
grows a multi-fault entry is the other half of this principle failing.

BP5, BP5-dip90, BP7 and BP8's regression compsets (`test.bp5.qdc`,
`test.bp5.qdc.dip90`, `test.bp7.qdc`, `test.bp8.qdc`) all have `ntotft = 1`.
Three multi-fault bugs in this project's history were found by running a new
case by hand; none were caught by the gate, because nothing in it ever
exercised `ntotft > 1`. `case_input/bp1002.qdc.2500` is the two-fault
step-over that now carries this: it has a reference under
`reference/bp1002/cycle0/`, a row in `tests/e2e/cases.py`, and a physical
invariant in `tests/unit/test_physical_invariants.py` asserting the seed
lands on fault 0 and only fault 0.

It sits in the **full** tier, not the every-push fast tier, because the event
is 3821 steps and ~2600 s on 3 ranks. So every-push CI still has no
`ntotft > 1` case. That is a known, costed gap, recorded in the CASES table
beside the row — not an oversight.

**How to apply**: keep a multi-fault case in the gate with a regenerable
reference and at least one assertion that distinguishes the faults from each
other, so a fault-aliased or misrouted read fails rather than merely moving a
number.

**Tier**: mechanical —
`python3 -m pytest tests/unit/test_physical_invariants.py -k multifault`.
The earlier citation here pointed at
`tests/contract/test_gate_set_has_multifault.py`, deleted in 2ccafe4; a rule
whose stated check does not exist is unenforceable, which is the failure mode
rule 12 is itself about.

---

## 8. Reference results are read-only, only grow, and every file has a reader

**Instantiates**: G8.

`reference/<benchmark>/*` are oracles (rule 3):

- Never edit or regenerate an existing gold file as a side effect of
  unrelated work.
- A new gold file is only committed after a run whose comparison verdict you
  have personally read — the first reference for a new benchmark is, by
  definition, unverified against any prior oracle. Say so explicitly in the
  commit: "first reference for X, not yet independently verified" rather
  than implying it carries the same evidentiary weight as an
  established zero-diff history.

**8a. Every reference artifact needs a test that reads it.** Gold that
nothing asserts on is dead weight masquerading as a safety net — an oracle
nobody consults would not have caught the incidents rule 3 fixed. Before
committing a new file under `reference/<benchmark>/`, confirm something
in `tests/`, `scripts/`, `tests/e2e/test_benchmarks.py` or `tests/e2e/cases.py` actually names
it (a literal string, or a dynamically constructed one — see
`tests/contract/test_reference_gold_is_referenced.py`, which checks both).

**Tier**: mechanical for read-only-and-growing (checksum/mtime diff on
`reference/` before vs. after a change should show only additions in the
change's own benchmark subtree); mechanical for 8a —
`python3 -m pytest tests/contract/test_reference_gold_is_referenced.py`.

---

## 9. Build/environment claims must be verified on the target host, not assumed

**Instantiates**: G9.

Known landmines, do not rediscover them the hard way:

- `src/makefile` hardcodes `FC = mpif90.openmpi` and `SCALAPACK_LIB = -L${LIB}
  -lscalapack-openmpi` for several targets. These names are Debian/Ubuntu
  `update-alternatives` symlinks; on a host without that alternatives system
  (or with a different MPI installed), only bare `mpif90`/`mpirun` exist.
  Before assuming a new host works, run `which mpif90.openmpi mpirun.openmpi`.
- **Build with `EQQUASIROOT=<root> MACHINE=utig make -C src`, then `mv
  src/eqquasi bin/`.** Bare `make` fails by design: `src/makefile` has an
  explicit `$(error MACHINE is not set. ...)` when `MACHINE` is unset — this
  is G2's "fail loudly" principle done correctly, not a bug to route
  anywhere. Do not silently export a default `MACHINE` to make a build
  succeed on a host it was never verified on.
- `local` target's `SCALAPACK_LIB` points at a system-wide
  `-lscalapack-openmpi`, which does not exist if MUMPS was built locally with
  its own bundled ScaLAPACK — check `./mumps/build/local/lib` first.
- `README.md`'s dependency list omitting `netCDF4` (while `scripts/case.setup`
  does `import netCDF4 as nc`) is a doc-vs-code drift, not a rule violation to
  fix here — route to `sophia-okafor`.

**Tier**: mechanical for `which mpif90.openmpi` / `ldconfig -p | grep
scalapack-openmpi`; judgment for whether a given host's fix belongs in the
shared makefile vs. a local override.

---

## 10. Docs move with the code, in the same change

**Instantiates**: G10.

When a compset, script, or file is added, moved, or renamed, `README.md`'s
matching section (e.g. "Currently supported compset includes") must change in
the same commit, and any pre-existing typo one line above your own edit
should be fixed while you're there.

**Tier**: mechanical for a cited typo (`grep -n <typo> README.md`); judgment
for "is the compset/doc list still complete."

---

## 11. Compare like with like: pin the version, vary one thing, stress the assumption

**Instantiates**: G11.

Three ways a comparison between two runs has produced a wrong conclusion here:

- **Different binaries.** A BP8 "resolution regression" was chased for hours
  before `runInfo.json` showed one run was v1.4.5 and the other v1.4.8 — the
  older predated a pore-pressure cell-area fix. `runInfo.json`'s `version`
  field (`src/library_output.f90:540`, written by every run — see
  `tests/contract/test_repo_hygiene.py::test_version_has_a_single_source_of_truth`
  and `::test_declared_version_matches_a_git_tag`) exists exactly to make this
  checkable: **read both runs' `runInfo.json` before trusting a diff between
  them, and check for a quantity that physically cannot depend on the varied
  parameter** — an unexpected change there is the tell that something besides
  the intended variable moved.
- **More than one variable changed.** The same investigation compared dx=50
  against dx=10 while the fault-normal extent and mesh grading also differed,
  and concluded refinement made things worse. It hadn't. Before citing a
  result as "X caused Y," diff the two `user_defined_params.py` (or
  `runInfo.json`) and confirm only the one intended parameter differs.
- **Verifying at a configuration that happens to satisfy an assumption,
  rather than stress it.** Multi-fault meshing was verified at dx=1000 with a
  5 km fault offset — 5 cells exactly. At dx=2000 the offset fell between
  mesh lines and one fault silently got zero nodes (rule 2a). A
  configuration where an integer-multiple assumption is exact is the one
  configuration least likely to expose a bug in that assumption; pick one
  that forces a remainder. Rule 12 is the same argument applied to the gate
  itself: the gate's own case selection can satisfy-not-stress just as an ad
  hoc comparison can.

**Route**: a `runInfo.json`-detectable mismatch is grounds to stop and reopen
the comparison, not to route anywhere — it means the comparison hasn't
happened yet.

**Tier**: judgment throughout — the *inputs* (`runInfo.json`'s version field)
are already mechanically produced and checked for internal consistency by
`test_repo_hygiene.py`, but *using* them to validate a specific comparison
between two runs is a discipline, not a static repo invariant, and cannot be
generically automated.

---

## 13. Known Fortran landmines in this codebase must not recur

**Instantiates**: G12. (13a also instantiates G6 — see the cross-reference in
rule 6.)

Three conventions this project adopted the hard way, each with a permanent
regression guard at `tests/regression/test_code_convention_landmines.py`:

- **13a. No `if (ntotft == 1) ... else ...` branching.** The fault-node
  engine was made neutral to `ntotft` in git `e5364d1`/`b6010e2` precisely so
  that `ntotft = 1` (every current oracle) exercises the identical code path
  `ntotft > 1` does. A branch here silently un-does that and leaves the
  multi-fault path unexercised again — the same failure mode rule 12 guards
  from the gate side. One code path, `ntotft = 1` the degenerate case.
- **13b. `nint` is not the intrinsic.** `globalvar.f90` declares `integer,
  parameter :: nint = 8`, shadowing the intrinsic function. Calling `nint(x)`
  elsewhere gives "Unclassifiable statement", not a type error — use
  `int(x + 0.5d0)`.
- **13c. Never open a data file on unit 6.** Unit 6 is stdout; opening a
  file there interleaves or clobbers every `write(*,*)` diagnostic.

**Tier**: mechanical, all three —
`python3 -m pytest tests/regression/test_code_convention_landmines.py`. All
three pass today (verified 2026-08-12); the guard is against regression, not
an open defect.

---

## 14. Shared-checkout discipline: worktrees, not branches; explicit paths, not `-A`

**Instantiates**: G13.

This repo runs several agents in parallel against one checkout and one
compute host.

- **Work in your own `git worktree`, never `git checkout` a branch in the
  shared checkout.** Branching in place has put two commits on another
  agent's branch here before. `git worktree add ../seas_bp10_eqquasi.<topic>
  <branch>` costs one command and cannot collide with anyone else's checkout
  state.
- **Stage explicit paths, never `git add -A` / `git add .`.** `git add -A`
  has swept another agent's unverified work-in-progress into a commit that
  wasn't about it. `git status` before staging, then `git add <path>
  <path>...` for exactly the files your change touches.

**Route**: a commit found to contain unrelated staged changes → the agent who
made the commit, to split it; if that's not tractable, `victor-reyes` for a
broader audit of what else landed unintentionally.

**Tier**: judgment — both are about what an agent chooses to run, not a
static property of the tree, so neither is mechanically preventable from
inside the repo. (A `pre-commit` hook rejecting `git commit` when the branch
under `.git/HEAD` differs from the worktree the agent was assigned is a
theoretical mechanical version of the first; not implemented, and it would
need agent-identity information this repo has no way to obtain.)

---

## 15. Operational safety on the shared compute host

**Instantiates**: G14.

All work here runs on one shared 64-core box (`cotopaxi`, per
`reference/bp8/runInfo.json`). All simulation artifacts belong under
`work/` (gitignored, see `tests/contract/test_repo_hygiene.py`), and launches
should default to `mpirun -np 1`, `OMP_NUM_THREADS=1` (already pinned by
`scripts/case.setup`'s generated launchers — see
`test_repo_hygiene.py::test_generated_launchers_pin_openmp_threads`) unless a
specific run has been coordinated with whoever else is on the box.

Three shell footguns that have each cost a run here, not specific to this
project but specific to how badly they interact with a shared long-running
box:

- **`pkill -f <pattern>` can match the issuing shell itself and kill it**,
  taking a long run down with it. Kill by PID: find it with `readlink
  /proc/<pid>/cwd`, not by pattern.
- **`echo 1 > file` redirects `echo`'s own stdout, it does not pass `1` as an
  argument to whatever reads `file`.** Depending on quoting, this can write
  an empty or wrong value where a control file was expected. Use `echo "1" >
  file` and check the file's contents afterward, not just the exit code.
- **`./case.setup` regenerates `run.sh` from scratch on every run**,
  discarding any manual edits to it (including a hand-pointed `mpirun.openmpi
  -np N <path-to-binary>`). Re-apply launcher edits after every `case.setup`,
  or template them into `scripts/case.setup`'s `create_run_sh` instead of
  editing the generated file.

**Also**: an agent must not end its turn waiting on a background job — there
is no mechanism that wakes it back up. Either block synchronously on the
result before finishing, or make it explicit to whoever picks up the thread
next that a job is still running and where its output will land.

**Tier**: judgment throughout. `case.setup` regenerating `run.sh` is a
verifiable fact about the script (`grep -n "def create_run_sh"
scripts/case.setup`), but "don't `pkill -f`," "quote your `echo`," "don't
strand a background job" are behavioral norms this file cannot enforce from
inside the repo — pretending otherwise would be worse than leaving them as
guidance.

---

## 16. A subagent's or audit's finding is a hypothesis until checked against source

**Instantiates**: G15.

Two examples from this project: an audit reported a "critical" friction-law
drift where the spec's own equation reads `ln(V*θ/D_RS)` with a reference
velocity and the code was in fact correct; separately, an agent's
nucleation-length estimate was off by 33×. A finding — including one in this
file, or one produced by the process this file describes — is evidence to go
check, not a conclusion to act on. Verify against the cited source (the
paper, the equation, the file:line) before it drives a fix, a routed report,
or a rule change.

**Tier**: unenforceable by construction — this is a rule about how a reader
treats a claim, which nothing in the repo can inspect. No command makes this
checkable; the only mitigation is the habit itself. Listed here rather than
omitted because a rule book that only states what it can check would hide
exactly the failure mode this rule is about.

---

## 18. The whole workflow must work for every example

**Instantiates: G7** (a configuration doesn't exist to the safety net until it's
registered where the net looks) and **G3** (a gate is only real if it runs and
fails loud).

Create the case, run it, post-process it. All three steps must work for **every**
compset in `case_input/`, and above all for the runs frozen under `reference/`.
A gold reference nobody can regenerate or plot is a file, not a reference.

Concretely, for each example:

- `create.newcase <dir> <compset>` succeeds and `./case.setup` produces a
  runnable case;
- the solver runs it to its own exit criterion, not a step cap;
- **every post-processing utility runs on the result** and produces a figure,
  or says clearly why it does not (BP8 is aseismic, so an empty rupture-time
  plot is the correct answer -- that is a message, not a traceback).

*Incidents, 2026-08-12, all found by pointing a tool at a benchmark it had never
been run against.* `plotOnFaultVars` failed on BP8 because it read `global.dat`
with `np.loadtxt`, which chokes on the section 4.2 field-name line that BP5 does
not have -- four utilities had the same defect and now all route through
`seasio.read_array`. `plotStations.py` was written against the BP8 column layout
and mislabelled the 9-column BP5/BP7 one throughout: slip rate is **linear**
there, not log10, and column 7 is effective normal stress, not pore pressure.
`plotOnFaultInitals` had never run at all. `plotAccumulated` and `plotProfiles`
were broken by a missing helper and a dead `pdf2image` import.

None of these were caught by a test, because the tests exercised one benchmark
each. The cheap guard is to run the utilities across all of them: differences in
column count, file extension (`.txt` against `.dat`), header presence and cycle
layout are exactly what a single-benchmark test cannot see.

Related: rule 8a (every reference file has a reader) covers whether a gold file
is *read*; this rule covers whether the pipeline that produced it still *works*.

## 19. Drive runs through the workflow, never the binary directly

**Instantiates: G3** (a gate is only real if it runs at the right
configuration) and **G18** (the whole workflow must work for every example).

To change what a run does, change `user_defined_params.py`, then `./case.setup`,
then `bash run.sh`. Do not invoke `mpirun ... bin/eqquasi` by hand.

`run.sh` is not a convenience wrapper. It owns the cycle loop: it writes
`currentcycle.txt`, moves each cycle's output into `cycle<i>/`, copies the
restart files back to the case root between cycles, and invokes the
post-processing. Calling the solver directly silently skips all of it.

*Incident, 2026-08-12.* The BP5 full-cycle gold was produced by running the
binary directly. Its output therefore landed flat in the case root instead of in
`cycle0/`, so it did not look like any other run; `plotOnFaultVars` then failed
on it because `user_defined_params.py` was not where a cycle directory would
have had it, and the file had to be copied in by hand to make the gold
plottable at all. Every conclusion drawn from that run was still valid -- but
the artifact was shaped unlike anything the workflow produces, which is the
opposite of what a reference should be.

The same applies to the parameters: edit the compset, do not hand-edit the
generated `model.txt` or `run.sh`. `./case.setup` regenerates both, so an edit
to either is discarded the next time anyone runs it (see rule 15).

**One parameter scheme, not two.** Every compset inherits from
`scripts/defaultParameters.py`, and the default must cover every case the code
supports -- otherwise a compset that needs something the default lacks builds it
by hand, and there are two ways to express the same thing.

*Open example.* `defaultParameters` carries `ntotft` and `faultgeom`, so it is
multi-fault aware, but allocates `on_fault_vars` as a single-fault
`(nfz, nfx, 100)`. `test.stepover.qdc` and `bp1002.qdc.2500` therefore allocate
their own 4-D `(ntotft, nfzMax, nfxMax, 100)` and fill it themselves. Two
schemes for one thing. The fix is the same `ntotft`-neutrality already required
of the Fortran (rule 13): allocate 4-D always, with `ntotft = 1` the degenerate
case.

**Run what the workflow puts on PATH, and know which binary that is.**
`install.eqquasi.sh` prepends `bin/` to PATH, but a stale EQquasi installation
elsewhere on PATH will shadow it and fail confusingly -- on this machine
`which -a eqquasi` finds a second copy under a different home directory that is
missing `libdmumps-5.4.so`, so `run.sh` dies with a shared-library error that
says nothing about the real cause. Check `which -a eqquasi` before concluding a
run is broken.

## 17. CI is the gate, not the local suite

**Instantiates: G3** (a gate is only real if it runs, fails loud, and is
evaluated at the right configuration) and **G15** (a finding is a hypothesis
until checked against the source).

A red CI run stops work. Do not commit, tag or release on top of one. Check it
after pushing -- `gh run list --repo dunyuliu/EQquasi --limit 5`, and
`gh run view <id> --log-failed` when it is red.

*Incident, 2026-08-12.* CI went red on the `mfault` merge at 03:41 and stayed
red through four tagged releases (v1.5.0, v1.6.0, v1.7.0, v1.7.1, v1.7.2)
because the local suite was green and nobody looked. The failure was real: the
merge added a `nid_fault` dimension to every on-fault netCDF write, so runs
produced `(1, nz, nx)` against gold's `(nz, nx)`, and `tests/e2e/test_benchmarks.py` compares
with xarray's `identical()`, which rejects differing dimensions.

**A comparator must not normalise away the difference it exists to detect.**
`scripts/plotAgainstGold.py` reported "matches gold" throughout, because its
`load_field` keys on `(variable, fault index)` and maps a 2-D and a 3-D array
onto the same key -- the shape difference was consumed before any comparison.
Two comparators disagreed for fourteen hours and the more permissive one was
believed.

Enforced in part by `tests/contract/test_gold_netcdf_shape.py`, which asserts
every gold fault snapshot carries `nid_fault` and stores its variables 3-D --
the cheap version of the CI signal, failing in seconds rather than after an
eleven-minute pipeline. `plotAgainstGold.py` now compares stored dimensions
before reshaping and refuses outright when they differ.

## Unenforceable rules (as written)

| Rule | Why it can't be checked today | What would fix it |
|---|---|---|
| 3 (bit-exact tolerance) | RESOLVED 2026-08-13. `tests/e2e/cases.py` compares at `rtol = 1e-9` for every benchmark, netCDF included. Bit-exactness is deliberately NOT the bar: two runs of the same case, same binary, same host differ by ~2e-14 from MPI reduction ordering alone, so `max\|diff\| = 0` fails on noise it cannot distinguish from a regression | — |
| 4 (model.txt positional contract) | No schema file pairs `case.setup` writes with `read_input.f90` reads — verifying order-correctness means manually diffing two files | A generated/shared schema (YAML or a Python list of `(name, type)`) consumed by both a `case.setup` codegen step and a Fortran read-order check |
| 6 (additive-only bp branches) | No script maps byte ranges of `if (bp == N)` blocks to a diff and rejects edits inside them | An AST/regex-based check scoped to `src/faulting.f90`, `src/library_output.f90`, `src/solveTimeLoopMUMPS.f90` |
| 7 (compsets.txt registration) | `compsets.txt` is not read by any script — it can drift from `case_input/` silently | `diff` the compset directory listing against `compsets.txt` in CI |
| 11 (comparison discipline) | "did you actually check `runInfo.json`, vary one thing, and pick a stressing configuration" is a property of an investigation, not of the tree at any commit | None generically; a per-investigation checklist is the closest available mitigation |
| 14 (worktree/`git add -A` discipline) | Which worktree an agent is "assigned" to and what it intended to stage are not repo state | A `pre-commit` hook could reject a commit whose staged-file set is a superset of `git diff --name-only HEAD` for files with a very recent unrelated mtime, but this is heuristic, not exact |
| 15 (operational shell footguns, background jobs) | These are properties of a command an agent chose to run, not of the repository | None from inside the repo; the mitigation is the norm itself |
| 16 (subagent/audit finding verification) | Whether a reader checked a claim against source is not observable from the repo | None — norm only |

---

## Rule-book health

**2026-08-12, structural pass.** The coordinator's original brief claimed no
rule book existed; that was a tooling error on their end (a glob that errored
silently read as "nothing found"), corrected after this file had already been
updated in place — which was the right call regardless of why it was needed.

That same day, at the owner's direction, this file was restructured
general-first/specific-second (Part I / Part II above) so it teaches the
transferable principle before it warns about the local instance, instead of
reading as a flat list of gotchas. No rule was culled, no incident was
dropped, and no rule number changed — rules 7 and 12 now sit adjacent in the
document because they share principle G7, which is a change in physical
order only; their identities (`#7`, `#12`) are exactly what they were before,
including in the three test files that cite them by number.

Earlier the same day: this file was found already seeded (dated 2026-08-02,
describing BP8 and the tiered `tests/` suite as future work) — not "no rule
book yet." The project has since reached v1.6.0: `test.reference.results/`
was replaced by `reference/<benchmark>/`, the
`tests/unit`/`contract`/`regression`/`e2e` suite now exists alongside
`pytest -m e2e`/`tests/e2e/test_benchmarks.py`, the fault-node engine was made
`ntotft`-neutral, and two of rule 2/3's three originally-cited violations
were fixed and are now mechanically guarded by
`tests/contract/test_release_gate.py`. That pass updated the stale paths and
violation claims in place, added rules 11-16 for incidents not previously
covered, and added three new checks (rule 8a, 12, 13) under `tests/contract/`
and `tests/regression/` rather than inventing a parallel mechanism. No second
rule book or stub was found anywhere in the tree, on either pass.

**Gap carried forward, not closed**: this file has no living status board
(starter invariant "a living status board, re-checked on a schedule"). The
closest thing today is `tests/README.md`'s tier descriptions and this file's
own per-rule "verified <date>" notes, which are not re-checked on any
schedule and are not a single place recording every open item. Not built in
either pass because it was out of the scope given (seed/enforce/audit/
reorganize the rule book, not create a new tracked artifact) — flagged so
it isn't silently forgotten.

Starter-set rules not carried over at original seed time: rule 9 (cheap
targeted check before expensive run) and the provenance/benchmark starter
rule were folded into rule 3 rather than kept generic, since this repo's
expensive-run hazard (regenerating gold references) already has a named gate.
