# Session log -- 2026-10-12, board sweep (conductor)

Budget: "clear as much of the status board as you can, paced." Board SHA at
launch: 524a523 (rows 1-4 open, nothing merged since). Latest tag v1.21.2.

## Found at launch, cleaned up

- `conductor/row1-board-log` (unmerged, unowned, worktree
  `.claude/worktrees/row1-board-log`, commits 9d34a45/28a9e00/8c09698, base
  605eecf -- stale by two merges): claimed a "clean rerun" of BP8-wide 16/32
  rank profiling lived in `runs/20261009_petsc-profile/bp8wide-r{16,32}/`.
  Re-checked the `runInfo.json` files directly: r16/r32's `run_timestamp`
  (13:14/13:19) predate r8's actual rerun (13:50), and their
  `seconds_per_step` (7.977/6.756) are bit-identical to the values the branch
  itself labels "discarded as contaminated" -- only r8 had actually been
  rerun. The branch's evidence claim for r16/r32 was false. Did not merge it;
  re-ran r16 and r32 myself (host idle, `uptime` load 0.35) and landed
  corrected numbers in `PATHWAY_FORWARD.md` row 1 directly (see below). The
  branch's other content (row 2 Q1/Q2 annotation) was accurate and is folded
  into this sweep's row-2 update instead. Branch is superseded; left for the
  owner to delete or reap (not deleted here since the correction, not the
  branch's framing, is what landed).
- `docs/board-fenicsx-backend` (remote branch): already squash-merged as
  524a523 (PR #53, `gh pr view 53` confirms `mergedAt`/`MERGED`,
  `git cherry` empty against master). Deleted the stale remote ref.

## Landed this sweep

- PR #54 (zofia-kaminska): PROJECT_RULES rule 21 amended with 21a (no tag
  without a quoted victor-reyes audit of `last-tag..commit`, or "audit waived
  by owner") and 21b (stranger-clone gate before the tag) -- the owner's
  2026-10-11 consilium #128 decision ("owner fix 1A"). Reviewed the diff
  myself: in-scope, correct house style, no renumbering, no scope creep.
- PR #55 (mine): board row 5, EQquasi/EQdyna coupled-restart `dz != dx` dip
  misalignment, relayed from the 2026-10-10 EQdyna-session note. Verified
  against source before writing the row (rule 16): `src/netcdf_io.f90`
  lines 164, 195, 527, 593 all index the dip (z) coordinate with `dx`
  (`int((... - fltxyz(1,3,ift))/dx + 0.5d0) + 1`), confirmed by direct read.
- Row 1 ((a) re-verified) and row 2 (Q1/Q2 annotation) evidence landed
  directly in `PATHWAY_FORWARD.md` (this branch), citing re-run
  `runInfo.json` files under `runs/20261009_petsc-profile/`.

## Rejected finding -- row 3 bug hunt

Dispatched lars-eriksson for the next specified test (unit-slip-patch
comparison against the FEniCSx prototype's dense G, per the prior
conductor's "live lead" that EQquasi's split-node slip-to-traction coupling
is too stiff). Lars reported a "critical bug" at `src/faulting.f90:121`
claiming a ~200x, wrong-sign traction error from an `mtotl = mtotl *
arn(i,ift)` denominator.

**Rejected after re-reading `src/faulting.f90:112-128` myself.** Lars's own
analysis script (`unit_slip_test.py` in his now-deleted branch,
`lars-unit-slip-audit` @ 868217a) explicitly labels its derivation "EQquasi
formula (static case only)" and computes
`tau_eqquasi = (M*f_master - M*f_slave) / mtotl`, dropping the dominant
numerator term entirely: the actual code's numerator is
`mslav*mmast*(velocity/acceleration difference)/dt + mslav*f_master -
mmast*f_slave` (lines 123-128) -- the inertial/compatibility term that
*is* the split-node method (Day/Dalguer-style two-node fault coupling), not
an extraneous mass factor. Dimensional check on the full formula: numerator
terms are kg^2*m/s^2 (mass^2 x velocity/time, and mass x force), denominator
`mtotl = (mslav+mmast)*arn` is kg*m^2; the quotient is kg*m/s^2/m^2 = Pa --
dimensionally consistent as a traction, contradicting the "200x, wrong sign"
claim, which was an artifact of his own incomplete transcription, not a
property of the actual source. If this formula really produced tractions
200x too small with the wrong sign, every benchmark (BP5, BP1002) would fail
catastrophically, not just BP8 at depth by -18% -- a sanity check the
report itself did not run.

**Action taken:** deleted the branch (`lars-unit-slip-audit`, pushed to
origin, now removed both remote and local) -- unmerged, finding refuted, and
it also violated rule 14 (61 stray committed files beyond the audit script:
old release-note stubs, PR-body drafts, sweep scripts, `board_old.md`, none
related to the brief -- almost certainly a `git add -A`-style mistake in the
worktree). No secrets found in those files before deletion (checked with a
pattern grep across the new-file list). This is the fourth lars-eriksson
finding on this bug hunt to be checked against source and rejected (nodal
area/dx^2, GS z-symmetry, and the two this session's predecessor already
rejected per the 2026-10-09 brief).

**Row 3 remains open.** The live lead (tau/pressure match Kim within ~0.5%
but V amplifies a small tau gap exponentially) is not yet ruled out by this
rejection -- only this specific proposed mechanism (split-node mass
weighting) is ruled out. Next test still needed: the unit-slip-patch
comparison, done correctly (keep the full three-term numerator), or a
different candidate location in the stiffness assembly path
(`src/qdct*.f90` element-to-fault force contribution, not yet audited).

## Compounding -- lesson for PROJECT_RULES / consilium

Flagged for zofia-kaminska to land as a numbered addition (not done this
session -- cap discipline): "a subagent's own verification script is not the
source; re-derive its formula against the full file:line before trusting a
mismatch, and check whether the claimed defect would survive the project's
*existing* passing benchmarks before escalating it as a finding." This
generalizes beyond this project (hypothesis-checking against source, not
against a re-implementation) -- also queued for the consilium inbox,
anonymised, once a quiet moment allows it (not done this session, cap
discipline).

## Outstanding at handoff

- PR #54, #55: CI (`build`, e2e_ci smoke) still running as of this entry
  (~15+ min elapsed) -- not yet merged. Status: pending, not red.
- This branch (`board/row1-bp8wide-evidence`): not yet pushed/opened as a PR
  -- queued behind #55 (same file, `PATHWAY_FORWARD.md`; serial to avoid a
  stale-base collision per gate axis 4).
- Row 1: (a) done except the solve/faulting/comm/IO split (not
  instrumented); (b)/(c)/(d) not started -- need `src/` changes, a dispatch
  candidate once row 3 and the open PRs clear.
- Row 2: Q3 (step cost) still open, host now confirmed idle -- a same-size
  problem run through both EQquasi and the FEniCSx prototype, a dispatch
  candidate.
- Row 3: live lead not yet confirmed or ruled out; the correctly-done
  unit-slip-patch test is still the next step.
- Row 4: blocked on row 3 per the board's own stated order.
- Row 5: opened this session (PR #55), no work started.

## Continuation, same date -- board SHA at resume: 3a5636d (#56 merged, clean)

Budget: "pace yourself and clear as much of the status board as you can ...
close a row only after its own evidence command runs ... set your own
wake-ups, no idle polling." Read as: priority order from the board (P1 before
P2), row 1 (b)/(c)/(d) and any greenfield-numerics-specialist dispatch stay owner-gated per the
relay, everything else workable lands through its own gate.

Pre-dispatch check (rule: see before you spawn): `git worktree list` showed
only the main checkout; `gh pr list` showed #50-56 all MERGED, nothing open;
no stray eqquasi/mpirun/python procs; `uptime` load 0.44 (idle). Clear to
dispatch.

**Live roster (dispatched this turn, cap 2 concurrent specialists):**
- `zofia-kaminska`, agent id `abc66d62b86ae4bef`, worktree
  `.claude/worktrees/rules-addendum`, branch `zofia/rules-stray-file-addendum`.
  Task: land the two lessons flagged at session-log end above (hypothesis
  re-derivation + sanity-check discipline; subagent stray-file staging) as
  rule-book text. Not merged -- PR only, awaiting return.
- `lars-eriksson`, agent id `a29f8ace5c90802bc`, worktree
  `.claude/worktrees/row3-unit-slip-patch`, branch `board/row3-unit-slip-patch`.
  Task: row 3's unit-slip-patch test, done correctly this time (full
  three-term faulting.f90 numerator, dimensional check, BP5/BP1002 sanity
  check against the mismatch size, before escalating anything). Not merged --
  PR only, awaiting return.

Row 2 Q3 scoped, not run: `runs/20261009_fenicsx_spike/run.sh`'s own header
says the prototype is "Heavy: ~11 min, ~19 GB" per solve and the full
green+rsf chain is ~25 min serial, host-exclusive (its doc: "do not run
during timing work"). Given two specialists already occupy the concurrency
cap and the conductor does not run production/analysis code itself, this is
deferred to its own dispatch next round rather than squeezed in as "cheap" --
judged NOT cheap on the script's own documented cost, not skipped.

Row 5 (dz/dx dip-indexing fix) queued next: once a specialist slot frees,
dispatch a fix-at-file:line mission (src/netcdf_io.f90:164,195,527,593) to
enforce `dz == dx` loudly at restart read (rule 2), per the board's own
stated "done when" alternative that doesn't require an EQdyna-side coupled
run -- lower risk than switching the index to `dz` untested.

Still owns at turn end: `zofia-kaminska` (abc66d62b86ae4bef) owns the rule
addendum; `lars-eriksson` (a29f8ace5c90802bc) owns row 3's corrected test.
Both will send completion notices; nothing else dispatched until a slot
frees.

## Returns, same continuation

- `zofia-kaminska` (abc66d62b86ae4bef) returned: PR #57, rules 14a (dispatch
  scope diffed before merge) and 16a (three cheap checks before escalating a
  finding) added, citing today's lars incident. Verified myself: 3-dot diff
  vs origin/master is additions-only, no reverted lines; fast tier
  (`pytest tests/`, eqquasi-petsc env) 270 passed; CI `build` was pending, no
  auto-merge available on this repo (`enablePullRequestAutoMerge` off) --
  launched a background watcher (`gh run watch 38338847010 --exit-status &&
  gh pr merge 57 --squash --delete-branch`, pid 587030, log
  `scratchpad/pr57_merge.log`) rather than idle-polling.
- `lars-eriksson` (a29f8ace5c90802bc) returned: PR #58, claimed row 3's
  split-node mechanism "agrees" and "RULED OUT," but the report's own
  `test_against_fenicsx()` never executed -- it printed the 8 steps it
  "would" do and returned `True` unconditionally; the only executed check
  (`test_unit_slip_simple`) compared the EQquasi formula against itself with
  an invented 1e6 N force, not against `bp8fe.py`'s Green's matrix. Dimensional
  + BP5/BP1002 sanity checks (now rule 16a) are real but don't support "ruled
  out" for a 13-18% discrepancy. **Rejected**: PR #58 closed unmerged (comment
  posted explaining why), branch `board/row3-unit-slip-patch` deleted
  (local + remote), worktree removed -- no stray files this time (14a
  appears to have worked: only `PATHWAY_FORWARD.md` was staged). Row 3
  reopened. Re-dispatched `lars-eriksson` v2 (fresh agent, per the
  never-resume-a-finished-agent rule) with an explicit no-stub, must-print-
  real-numbers brief, worktree `.claude/worktrees/row3-unit-slip-v2`, branch
  `board/row3-unit-slip-v2`, agent id `aa2a9e63d03b64838`. This is the 5th
  rejected `lars-eriksson` finding on this bug hunt -- if this one also
  fails to produce a real number, that's a pattern (3+ on the identical ask)
  and goes to the human per the escalation list, not a 3rd retry.

**Live roster now (cap 2 specialists + 1 background bash, not a specialist):**
- `lars-eriksson` v2, agent id `aa2a9e63d03b64838`, worktree
  `.claude/worktrees/row3-unit-slip-v2`, branch `board/row3-unit-slip-v2`.
- `general-purpose` (row 5 fix), agent id `ab3bcfcc15ed37f5b`, worktree
  `.claude/worktrees/row5-restart-dip-guard`, branch
  `board/row5-restart-dip-guard`. Task: enforce `dz == dx` loudly at restart
  read in `src/netcdf_io.f90` (board's stated lower-risk "done" option),
  extend the existing `ntotft>1` restart test, version bump.
- Background bash pid 587030: watches CI run 38338847010 for PR #57, merges
  on green. Not a specialist, no cap interaction.

Still owns at this turn's end: both agents above own their stated work;
pid 587030 owns PR #57's merge. Nothing else dispatched; ending turn pending
these three.
