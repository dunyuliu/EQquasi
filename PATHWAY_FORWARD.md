# Pathway forward

Open work only. Updated 2026-10-08, at v1.21.2.

Close an item by deleting it, and record what was learned where it will be
read again: `PROJECT_RULES.md` for a rule, `data/<bench>/README.md` for a
benchmark finding, a commit message for a fix. Closed rows and the findings
behind them (BP8 domain and boundary study, step-over science, caps, known
limits) are in `docs/dev/board_history.md`. Rules are in `PROJECT_RULES.md`;
build, run and traps are in `CLAUDE.md`.

---

## Open

1. [ ] **BP8-PW reference.** Decide it against the CRESCENT comparison; both
   10 m entries (EQquasi 1.21.0) were uploaded 2026-10-07. Check first: an
   immediate same-binary rerun differed by up to ~8 % at the 4 on-axis
   stations (same ranks, `OMP_NUM_THREADS=1`, `src/porepressure.f90`).
2. [ ] **BP8 convergence statement (optional).** On the wide box the residual
   vs Kim shrinks with cell size (along depth: -18 % at 50 m, -13 % at 10 m);
   domain size is converged (< 1 %). A 25 m wide-box run (~2 h, 8 ranks)
   would give a third point.
3. [ ] **Owner: two local branches with unmerged commits.**
   `docs/codify-pr-workflow` (2026-09-23, also on GitHub; a draft of the PR
   rule, likely superseded by rule 17) and `kai/bp8-simplify` (2026-08-07,
   local only; test tiering, likely superseded). Keep or delete.
4. [ ] **Owner: enable GitHub Pages** (Settings -> Pages -> Source: GitHub
   Actions) so `docs/user/` publishes on the next tag.
5. [ ] **Pre-existing tool bugs** (found writing `docs/user/`):
   `scripts/case.submit` and the `batch.hpc` it submits use the pre-1.13 case
   layout (the solver would stop with code 6); `scripts/plotAgainstGold.py`'s
   usage lists choices it does not accept and looks for a missing reference
   file for `test.bp5.qdc.2000`.
