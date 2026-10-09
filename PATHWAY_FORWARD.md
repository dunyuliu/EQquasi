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

1. [ ] **Owner: two local branches with unmerged commits.**
   `docs/codify-pr-workflow` (2026-09-23, also on GitHub; a draft of the PR
   rule, likely superseded by rule 17) and `kai/bp8-simplify` (2026-08-07,
   local only; test tiering, likely superseded). Keep or delete.
2. [ ] **Owner: allow tag deploys to Pages.** Pages is on (Actions), but the
   `github-pages` environment allows only `master`, so the tag-triggered
   deploy is refused. Settings -> Environments -> github-pages -> add a Tag
   rule `v*`; then re-run the latest docs workflow run.
