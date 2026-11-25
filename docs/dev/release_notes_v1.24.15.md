# EQquasi v1.24.15 — 2026-11-25

Tag target: the squash-merge commit of the release-notes PR on top of `b997272342aba877f1eeb61c8eb59cbe689f62d3` ("Audit minors: dip range, visible empty category, restart remedy copies first (1.24.15) (#288)"). `EQQUASI_VERSION = '1.24.15'` (`src/globalvar.f90:10`). Previous tag: v1.24.14 (`e11190e`). One version in range (rule 21c): `git log v1.24.14..origin/master -S'EQQUASI_VERSION' -- src/globalvar.f90` lists only #288.

Code diff `v1.24.14..origin/master` in `src/`, `scripts/`, `tests/` is exactly #288 (5 files, 37 insertions, 12 deletions: `scripts/generateFaultInterface`, `src/globalvar.f90`, `src/netcdf_io.f90`, `tests/e2e/test_benchmarks.py`, `tests/unit/test_fault_interface_dip.py`); the other commits in the range are board files.

## Behaviour a user can see

- **`generateFaultInterface` refuses a dip outside 0 < dip <= 90 (#288).** The guard now runs before the offset `dx/tan(dip)` is computed (`scripts/generateFaultInterface:101-106`), so a negative dip, a dip above 90 and dip 0 all stop with `par.dip=<value> is not a dipping fault`. Before, only `dip % 180 == 0` was refused and a negative dip silently built the same plane as its positive value. Only `bp1001.fdc.250` and `test.bp5.qdc.dip90.2000` reach the guard among shipped compsets, both with dip 90, so no shipped case changes. Shown by `tests/unit/test_fault_interface_dip.py`.
- **A restart missing `dip_node_spacing` now says to copy first (#288).** The message in `src/netcdf_io.f90:776-783` prints `cp fault.r.nc fault.r.nc.bak` before the stamping command; it said "a copy" while the command edited the original in place.
- **A deleted e2e category now fails by name (#288).** `category_params` in `tests/e2e/test_benchmarks.py:95-105` keeps one param `<category>-no-reference-holds-it` for a category no reference holds, and `test_category_matches_reference` asserts `no reference holds any <category> files`, so the category cannot vanish at collection. Collection on today's data is unchanged: 77 entries, 9 categories.

## Gate evidence

- Full tier (`tests/release_gate.sh`) on `7db7a95` (PR 288's final head; its `src/`, `scripts/`, `tests/` are identical to `b997272`'s): `runs/20261125_c45_gate_7db7a95/status.txt` reads `done PASS` (`clone rc=0`, `build rc=0`, `fast rc=0`, `e2e_full rc=0`, machine `conda-linux`). `e2e_full.log`: `145 passed, 11 skipped, 417 deselected, 3 warnings in 14500.08s (4:01:40)`. Fast tier in that gate: 417 passed, 0 skipped.
- All 11 skips are one reason: `SKIPPED [11] tests/e2e/test_benchmarks.py:127: no reference at .../data/test.bp1002.qdc.zone.2500/cycle0` (the open owner question on the board; the count matches the board's recorded 11, rule 3a).
- CI on `7db7a95`: run 45637120122, success.
- **The full-tier gate was not rerun on the tag SHA.** The tag goes on the commit that carries these notes; against the gated tree it differs only in board files and this file. On the tag SHA: CI (below) and the stranger clone.
- No performance claim is made by this release.

## Tag-SHA evidence

- CI on `b997272` (master head before the notes): run 45707597173, conclusion success (`build` success), https://github.com/dunyuliu/EQquasi/actions/runs/45707597173.
- CI on the notes-merge commit (the tag SHA) and the stranger clone (`runs/20261125_c45_stranger/`): recorded on the Release page and in the conductor's report, because a note cannot hold the run id of the commit that contains it.

## Audit

Code audit of `7db7a95` (PR 288's final diff; `git diff --stat v1.24.14..origin/master -- src scripts tests` is that diff alone, so it covers `v1.24.14..b997272` code), read-only independent code auditor, 2026-11-25: **0 blocker, 0 major, 3 minor.**

- Minor (a), `tests/e2e/test_benchmarks.py:104`: the fail-by-name fallback param carries the `e2e` mark but not `e2e_fast` or `e2e_ci`, so a deleted category shows only in the periodic full tier.
- Minor (b), `scripts/generateFaultInterface:119-121`: a leftover comment block still describes the old dip-0 check.
- Minor (c), `tests/unit/test_fault_interface_dip.py`: the new test adds `tests/e2e` to `sys.path` without undoing it; `monkeypatch.syspath_prepend` would.

Also reported: collection on today's data is identical (77 entries, 9 categories); no shipped compset path breaks (only `bp1001.fdc.250` and `test.bp5.qdc.dip90.2000` reach the guard, both dip 90).

## Audit findings and fixes

Applied here: none (all three minors are test or comment edits and ride in their own PR; none changes a behaviour a user sees).

Open follow-ups, not fixed in this release: minors (a), (b) and (c) above; the 11 reference-less bp1002 zone skips (owner question); issue 270 stays open as recorded in the v1.24.14 notes.

## Files added / removed

Added: this file. Removed: none.
