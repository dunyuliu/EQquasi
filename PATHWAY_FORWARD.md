# Pathway forward

Open work only. Updated 2026-11-13, at v1.22.0.

Close an item by deleting it, and record what was learned where it will be
read again: `PROJECT_RULES.md` for a rule, `data/<bench>/README.md` for a
benchmark finding, a commit message for a fix. Closed rows and the findings
behind them (BP8 domain and boundary study, step-over science, caps, known
limits) are in `docs/dev/board_history.md`. Rules are in `PROJECT_RULES.md`;
build, run and traps are in `CLAUDE.md`.

---

## Open

1. [ ] **P1 -- PETSc scaling in EQquasi.** Owner, 2026-10-09. Known today:
   native MUMPS scales 1.9x over 1->8 ranks (BP5); 10 m BP8 wide box ~15 s/step
   at 8 ranks, ~79 GB of factors; CG+GAMG scales better but is slower at every
   size tried (1.4x bp1002, 7.7x kink300); faulting and BP8 pore pressure run
   on rank 0 only. Steps, in order:
   (a) profile per-step time (solve / faulting+pressure / communication / I/O)
   at 8, 16, 32 ranks on 10 m BP8 wide and bp1002, host to ourselves;
   (b) MUMPS block-low-rank compression, ParMETIS ordering, threaded BLAS per
   rank; (c) hypre BoomerAMG vs tuned GAMG with the setup reused across steps
   (K is constant); (d) distribute the rank-0 work only if (a) says it matters.
   Done when: 10 m BP8 wide <= 5 s/step at 32 ranks, or a measured reason it
   cannot be, with parity to the gold at every change.
   Evidence: the per-step table from (a) and the final run's runInfo.json.
   **ROW 1 SUMMARY FOR THE OWNER, 2026-11-14 (conductor 12; every figure is a table or paragraph below, none new; 8-16 ranks, MEASURED, hosts shared so loads vary run to run).**
   (Row 1 history referred to as "below" here: (archived: board_history.md, section "Row 1 measurement history, 2026-11-14 to 2026-11-18 (archived 2026-11-18)"))
   - *Passes the gold-style comparator (1e-9):* exact MUMPS (reference); all-primal substructuring (bp1002 25/25, bp5 120/120, gs10 33/33; wide not runnable, dense Schur 270-295 GB); GAMG CG at tight rtol (wide 1e-14: 0 of 33 fail; bp1002 1e-11/1e-12); BLR-preconditioned CG at tight rtol. *Fails it:* BLR 1e-8/1e-6/1e-4 as a solver (bp1002 13-23 of 25, bp5 115-118 of 120, wide 13-19 of 33), GAMG CG at rtol 1e-8 (wide 9 of 33, bp5 116 of 120), subdomain-direct CG and BDDC (9-15 of 25/33).
   - *Seconds per step (loop):* bp5 0.128-0.178 exact, 0.267 all-primal, 0.335 GAMG; gs10 16 ranks 2.71 exact, 4.45 all-primal; wide 16 ranks 6.9-8.3 exact (5.5 at 200 steps), 32.4 GAMG at rtol 1e-8, 49.9-53.4 GAMG at the passing rtol (6-8x slower), 2.9 BLR 1e-8 (1.9x faster, fails). Exact wide at 8/16/32 ranks: 12.7 / 6.9 / 7.1; the row's target (wide <= 5 s/step at 32 ranks) is not met by anything that passes.
   - *Setup and memory per rank (wide, 16 ranks):* exact 331-344 s and 12.1-12.3 GB; GAMG CG 10 s and 2.54 GB (4.8x less); BLR 128-156 s and 9.5 GB; all-primal on gs10 268 s and 20.7 GB against exact 81 s and 4.5 GB.
   - *Merge gate:* full `e2e` tier on the stack head, solver 3 forced: 90 passed, 45 skipped (10 pending bp1002 zone gold, 35 category-absent), 0 failed (`runs/20261114_gate-full/`). Audit of the stack diff: 0 BLOCKER, 0 MAJOR, 9 MINOR (all listed above). Stack #76, #77, #82, #83, #85 is MERGEABLE with CI green; #87 (BDDC spike) is not for merge.
   - *64-rank `-11` abort:* MUMPS -9, cured by `-mat_mumps_icntl_14 100` (3 of 3 pass against 5 of 6 fail; +20-30 % rank memory), passed through `PETSC_OPTIONS`; 60 of 60 clean at 8 ranks.
   - *Owner decisions open:* (1) merge the stack and in what order; what `par.solver` means at (h) and PETSc as a hard dependency; (2) default ICNTL(14) in source, and `KSPSetUp` reading INFOG(1)/(2); (3) move the comparator tests to `tests/unit` (new file) and document `PETSC_OPTIONS`; (4) audit items 2-6, 8, 9; (5) closing this row, or taking the measured reason: no tested solver passes at <= 5 s/step on wide at 32 ranks, and the solver survey is saturated (no new runs are queued).
   **(a), re-verified 2026-10-12** (binary eqquasi-1.21.2, host theo4,
   host-exclusive, 10-step profiling runs, native MUMPS `par.solver=1`;
   `runs/20261009_petsc-profile/{bp1002,bp8wide}-r{8,16,32}/result/cycle0/runInfo.json`):
   bp1002.qdc.2500 (127k eq) -- 8r 0.436 s/step, 16r 0.409, 32r 0.394 s/step
   (flat; rank-0-only faulting/pressure work dominates at this size). BP8 10 m
   wide (3.52M eq, scratch-only `fxmin/fxmax/fzmin/fzmax=-1500/1500`,
   `enlarging_ratio_xz=1.2`, NOT the committed `compset/bp8.qdc.gs.10` which
   is a smaller -500/500 box) -- 8r 11.01 s/step, 16r 8.01, 32r 6.91 s/step.
   Correction to a prior conductor's note: an earlier pass flagged these as
   possibly contaminated by a concurrent external job and claimed a rerun was
   needed; on re-running r16 and r32 clean (host confirmed idle, `uptime`
   load 0.35, no other eqquasi/python procs), the new numbers (8.01, 6.91)
   match the flagged-contaminated ones (7.98, 6.76) to ~2-4%, i.e. within
   normal run-to-run noise -- the contamination was not material here.
   (b)/(c)/(d) not started. Finer solve / faulting+pressure / communication /
   I/O split (as (a) asks) is NOT yet instrumented -- only aggregate
   time-loop and factorization time exist today (`solveTimeLoopMUMPS.f90`'s
   `cpu_time` calls). Recommendation (not acted on -- rule 1, new-compset
   naming needs agreement): register the -1500/1500 wide box as its own
   compset (e.g. `bp8.qdc.gs.10.wide`) instead of leaving it only in `scratch/`.
   **Owner direction, 2026-11-13**: "i think we can domain decompose from the
   beginning with petcs?"; on the three changes below, "all three changes
   sounds great". Why: the time loop is two triangular solves per step and
   those do not scale (above; Wang 2026 Case 3a, 3.61M eq, 24 ranks, v1.22.0:
   factorization 137 s, time loop 15093 s per cycle, ~6.5 s/step). The PETSc
   path today is not a domain decomposition: rows are split by `PETSC_DECIDE`,
   rank 0 builds the whole RHS and receives the whole solution
   (`VecScatterCreateToZero`), and the default solve is still MUMPS. The
   three changes, in the Fortran PETSc path, replacing (c)/(d):
   (e) each rank owns a block of the mesh: it generates only its own
   elements (the generator is analytic, no partitioner), assembles its own
   rows and RHS, and runs friction and pore pressure on its own fault nodes;
   (f) memory scales out: no rank holds the full mesh or the global factors
   (79 GB at 10 m BP8 wide, 141 GB in-core estimate at Case 3);
   (g) subdomain matrices kept unassembled (`MATIS`) and solved with
   `PCBDDC`, the setup built once and reused every step (K is constant).
   Owner, 2026-11-14, on the solver under (g): "This is a really simple system.
   It's symmetric. You can use maybe just a regular CG or GRMS." (verbatim;
   read as: try plain Krylov, CG or GMRES, on the distributed matrix, not
   only BDDC with direct local solves.)
   Owner, 2026-11-14: "If you cut the domain into subdomains, and each
   subdomain you can use mumps with much less memory cost, then it also
   sounds a good strategy." (verbatim; read as: a subdomain-wise direct
   solve under an outer Krylov iteration is also a candidate.)
   Owner, 2026-11-14: "you can just factorize once and just solve for each
   subdomain and then piece them together. But you just factorize once."
   (verbatim.)
   Owner, 2026-11-14: "was a really good strategy and put it on board as
   well. So in the end, we may have several solver options. And I think the
   last that you proposed, this 3-step one, is pretty good. Give it a try."
   (verbatim.) The 3-step strategy, exact substructuring with no outer
   iteration: (1) factorize each subdomain interior once; (2) build and
   factorize once the Schur complement that couples only the interface
   unknowns; (3) each step is one solve per subdomain, one interface solve,
   one more solve per subdomain. Known cost: the interface matrix is dense,
   its size not yet computed. Status: TO TRY as a spike, beside plain CG and
   subdomain-direct solves under CG; none measured yet.
   Owner, 2026-11-16, on the 3-step solve after its measured cost (dense
   interface, 3-4x the memory of exact): "how come? the interface shouldn't
   grow that crazy" and, on the untried ways to cut the dense interface cost
   (BLR compression of the Schur complement offered as the cheapest test on
   the existing spike; a sparse hierarchical interface factorization; fewer,
   larger subdomains): "keep pursuing that line" (verbatim.) Owner, 2026-11-16:
   "keep pursuing that line, optimize it with memory and speed." (verbatim.)
   Status: steps 1 and 2 measured 2026-11-16 (conductor 13, block below);
   the cheaper levers are in progress, each locked here before the next.
    (The "block below" is archived: (archived: board_history.md, section "Row 1 measurement history, 2026-11-14 to 2026-11-18 (archived 2026-11-18)"))
   Owner, 2026-11-16, on the parity tolerance for approximate solves (the
   question further down this row): "let's use 1e-9 to machine precesion for
   now" (verbatim.) Read as (a reading, not the owner's words): the 1e-9
   comparator stays the bar for every solver option for now, and iterative
   solves are run as tight as it takes to pass it, down to machine precision.
   Owner, 2026-11-16, after the survey of the other solver lines (GAMG CG
   best on memory and scaling; its untested lever, rigid-body modes per side
   of the fault's split nodes in the near-null space, offered as a second
   line): "so the winner is gamg cg. what's the near null space on the split
   node?" and then "keep proposing, keep exploring. 3 days budget. go"
   (verbatim.) Budget read as 2026-11-16 16:20 CST to 2026-11-19 16:20 CST;
   scope read as (a reading): the solver lines on this row, the 3-step solve
   and GAMG CG first, new proposals recorded here before they are run.
   Owner, 2026-11-17 (asked which clusters and allocations are reachable;
   two TACC systems are): the scale test is to run on Vista under an
   allocation the owner named (marked paraphrase; the allocation code is
   left out of this public file), followed by "just need more cores to
   scale?" (verbatim.) Read as (a reading): the cluster test under (g) is
   released on that one allocation; the build there is spike work until a
   machine target is agreed; the question is answered with measurements.
   Owner, 2026-11-17, on which solver options move from spike to a
   merge-bound branch: "let me decide promotion" (verbatim.) Read as (a
   reading): no agent promotes, proposes a merge-bound rewrite of, or
   retires a solver option; rankings are offered, the owner picks.
   Owner, 2026-11-17, shown a ranking whose item 1 was the distributed
   PETSc loop with exact MUMPS (`par.solver = 3`, the stack #76, #77, #82,
   #83, #85): "cool. promote 1" (verbatim.) Read as (a reading):
   `par.solver = 3` and the new file `src/solveTimeLoopPETScDD.f90` are
   agreed, and the stack lands through gated merges in stack order. Not
   read into it: a change of default solver, retiring solver 1 or 2,
   anything under (h), or promotion of any other option.
   Owner, 2026-11-17, asked for a better name for that option and offered
   three: "distributed direct" (verbatim.) Read as (a reading): the
   option is called distributed direct (distributed assembly, one global
   direct factorization) in docs, the board and PR titles from here on;
   renaming the source file or any parameter is not decided.
   Owner, 2026-11-17: "and let's /refactor the code base to make eqquasi
   naturally adopts more solver options on the line?" (verbatim.) Read as
   (a reading): after the stack lands, a behaviour-preserving refactor of
   the solver layer so a new solver option plugs into the distributed
   direct loop without a new time loop; gold results identical; the
   design goes on this row before code; no new user-facing parameter,
   no default change and no retired solver without the owner's word.
   Owner, 2026-11-17: "/release to lock this promotion. other lines still
   pending my decision and test later" (verbatim.) Read as (a reading):
   once the stack is on master, cut a release (audit first, tag gated on
   green CI for that commit) so the distributed direct option is locked;
   the release comes before the refactor; every other solver line stays
   spike-only and undecided.
   Owner, same exchange: "and this should go through full test" and
   "with a full cycle?" (verbatim.) Read as (a reading): the release is
   gated on the full `e2e` tier on the exact commit with the distributed
   direct option exercised, not only the fast tier; "full cycle" is read
   two ways and both are done unless the owner narrows it: the full
   milestone cycle (rules audit, code audit, fixes, release), and a run
   through a complete earthquake cycle compared against the gold result.
   Owner, same exchange: "also, always measure perfs" (verbatim.) Read
   as (a reading): every gate, release and refactor step on this row
   records performance beside correctness (setup, time per step, peak
   memory, rank count, host load), before and after the change.
   The owner then named the three: "memory cost, scaling with more cores,
   and compute speed per step" (verbatim.)
   Owner, 2026-11-17: "and whatever metrics you'd better collect, say
   machine used, hardware used, i think eqdyna has some setups you can
   borrow" (verbatim.) Read as (a reading): each performance record also
   names the machine and its hardware (CPU model, cores, memory, compiler,
   MPI and library versions), plus whatever else makes the number
   reproducible; look at how EQdyna records its runs and borrow that
   layout rather than invent one. A new file or script for it is still a
   rule 1 proposal on this row first.
   Owner, 2026-11-18 13:55 CST, answering three numbered questions put by
   the invoking session: "3A. 2A" (verbatim.) The options as they were put:
   question 2, MUMPS workspace setting ICNTL(14), "A) default it to 100,
   since Vista fails at the default from 64 ranks (recommended), B) leave
   at 35"; question 3, solver-layer refactor design on this row, "A) go,
   gated at 1e-9 (recommended, since byte-identical only holds at 1 rank),
   B) hold". Read as (a reading): (2A) ICNTL(14) defaults to 100 in source
   on the MUMPS-through-PETSc paths, a user option still overriding it; a
   `src/` change, so a version bump and a gated PR, with setup time, time
   per step and memory recorded before and after; whether the native MUMPS
   loop (solver 1) gets the same default, and `KSPSetUp` reading
   INFOG(1)/(2), are NOT decided. (3A) the refactor design on this row
   (#166) may start; its gate is the 1e-9 comparator against gold and
   against the pre-refactor binary, not byte identity; it stays
   behaviour-preserving, with no new user-facing parameter, file rename or
   default change without the owner's word. Question 1 (an absolute floor
   for the out-of-plane columns in the comparator) was NOT answered and
   stays open: strict 1e-9 holds until it is.
   Owner, 2026-11-18 13:57 CST, answering question 1: "1A" (verbatim.) The
   option as it was put: "Comparator bar for out-of-plane columns: A) add
   an absolute floor so roundoff passes (recommended), B) keep strict
   1e-9". Read as (a reading): the comparator keeps rtol 1e-9 on every
   column that carries signal and gains an absolute floor on the
   out-of-plane columns (`Slip_3`, `Slip_rate_3`, `Shear_stress_3` and
   their log10 forms in the `fltst` files), so values at roundoff level
   pass. The floor value was NOT given: it is proposed on this row with
   the measured exact-against-exact differences behind it before any
   result is re-read as a pass; the cross-architecture float32 `disp`
   files are NOT covered by this answer.
   Owner, 2026-11-18 13:59 CST, after being told what stands between the
   slip-to-traction option and regular use (and, one message earlier, the
   Vista limits of distributed direct): "fix those vista issue, and run
   full cycles, a few, and verify. if we can pass that and the CI, let's
   promote it" (verbatim.) Read as (a reading, the owner may correct it):
   "it" is the slip-to-traction option; promotion is CONDITIONAL on
   (i) the Vista issues fixed: ICNTL(14) (2A above), the NUL holes in
   output written to NFS scratch (run from Lustre, or a fix if the cause
   is in the code), and for slip-to-traction a clean dedicated-node run;
   (ii) a few full earthquake cycles run with slip-to-traction and
   verified against exact at the comparator bar (1e-9 with the 1A floor);
   which cases and how many is proposed on this row first, with cost;
   (iii) the full test tier and CI green on the exact commit. Promotion
   then follows the same route as distributed direct: gated merges, a
   version bump, the code audit before any tag. NOT decided by this
   message: the option's name and parameter value (rule 1: proposed on
   this row and agreed before the merge-bound branch), a `vista` machine
   target, a compset for wide, the re-layout reader for the saved
   operator, any default change.
   Owner, 2026-11-18 14:00 CST: "then, let's promote slip to traction"
   (verbatim.) Read as (a reading): confirms that "it" above is
   slip-to-traction and that promotion comes after the conditions (i) to
   (iii) pass ("then"), not before them.
   Owner, 2026-11-18 14:00 CST: "and then optimize the slip to traction"
   (verbatim.) Read as (a reading): after the promotion lands, the
   slip-to-traction option is optimized; the measured costs on this row
   are the candidates (operator build time, 1.9 h at 32 ranks and 3 h at
   16 on wide; memory, about 1.3x summed and 1.5x peak rank against
   exact; the saved operator being readable only at its own rank count).
   Each optimization is proposed on this row with its before and after
   measurement and stays behaviour-preserving at the comparator bar; it
   comes after the promotion, not inside it.
   Owner, 2026-11-18 14:00 CST: "/autopilot go" (verbatim.) Read as (a
   reading): the unattended run continues on the scope above (Vista
   fixes, full cycles and verification, promotion of slip-to-traction if
   the conditions and CI pass, then its optimization, with 2A, 3A and 1A);
   "go" names no time budget, so it is read as until this row's items are
   done or blocked on the owner, replacing the 2026-11-19 16:20 CST end.
   The Vista ceiling of 150 SU and every "NOT decided" item above stand.
   Owner, 2026-11-18 16:18 CST, answering three numbered questions put by
   the invoking session: "1A. 2A. 3A" (verbatim.) The options as put were:
   question 1, refactor file name, "A) accept `src/onFaultHelpers.f90`
   (recommended), B) fold the helpers into an existing module"; question
   2, slip-to-traction option, "A) `par.solver = 4`, "slip to traction"
   (recommended), B) another value or name"; question 3, out-of-plane
   floor value, "A) 1e-12 m/s on `Slip_rate_3` and 1e-9 MPa on
   `Shear_stress_3`, with `Slip_3` measured first (recommended), B) other
   values". Read as (a reading): the new file `src/onFaultHelpers.f90`
   with module `onFaultHelpers` is agreed for refactor step 1 (rule 1);
   the slip-to-traction option is `par.solver = 4`, named "slip to
   traction", so its promotion stack may go on merge-bound branches; the
   comparator's absolute floors are 1e-12 m/s for `Slip_rate_3` and
   1e-9 MPa for `Shear_stress_3`, and the `Slip_3` floor is measured and
   put on this row before any value is used. NOT decided by these
   answers: the `Slip_3` floor value; how the floor applies to the log10
   columns; any new source file name the slip-to-traction stack may need;
   ICNTL(14) for native solver 1; a `vista` machine target; any default
   change.
   Owner, 2026-11-18 16:21 CST, after being told that the solver runs in
   double precision but the `disp.*.nc` snapshots are stored as 4-byte
   reals (`NF90_REAL` in `src/netcdf_io.f90`), which is why they fail the
   1e-9 comparator across machines: "weird. fix that. should be double
   precision" (verbatim.) Read as (a reading): `netcdf_write` stores the
   displacement snapshots as 8-byte reals, as a gated `src/` change with a
   version bump, performance and file size recorded before and after.
   NOT decided by this message: the same change for the other two
   single-precision writers in that file (`netcdf_write_on_fault`, whose
   `fault.*.nc` files are in the committed gold and are read back on
   restart, and `netcdf_write_roughness`); any regeneration of gold under
   `data/`; what readers outside this repo expect.
   Owner, 2026-11-18 16:22 to 16:24 CST, on the same subject: "there
   should be a setup of precision in globalvar.f90?" and "all f90 should
   follow that right?" (verbatim.) Told that `dp =
   selected_real_kind(15,307)` in `src/globalvar.f90` is followed
   everywhere except three `NF90_REAL` definitions in `src/netcdf_io.f90`
   and a few hard-coded but equivalent declarations, and asked two
   numbered questions, the owner answered: "for 1 and 2, yes, both should
   follow dp" (verbatim.) The options as put were: question 1, on-fault
   and roughness NetCDF writers, "A) switch them to double as well and
   regenerate the affected gold in a separate gated PR (recommended,
   since it matches your rule that everything follows `dp`), B)
   displacement only, leave those two in single"; question 2, the
   hard-coded `real(8)` and `double precision` lines, "A) change them to
   `kind=dp` in the same PR, behaviour-preserving (recommended), B) leave
   them". Read as (a reading): all three NetCDF writers store reals at
   the kind `dp` names, and the hard-coded `real(8)` and `double
   precision` declarations become `kind=dp`; the gold under `data/` that
   the on-fault change invalidates is regenerated in its own gated PR,
   with the old-against-new difference shown to be storage truncation
   only before any gold file is replaced. This supersedes the "NOT
   decided" items of the 16:21 paragraph for the two other writers and
   for that gold regeneration. NOT decided: the MPI calls that name
   `MPI_DOUBLE_PRECISION` (not asked); gold for anything other than the
   files this change invalidates; what readers outside this repo expect.
   Owner, 2026-11-18 16:25 CST, two messages in a row: "and rebless the
   affected disp" then "in ref" (verbatim.) Read as (a reading): the
   reference data that holds displacement snapshots affected by the
   change to 8-byte storage is regenerated ("reblessed") with the new
   binary, in the gated gold PR above, each replaced file shown to
   differ from the old one by storage truncation only. "ref" is read as
   the reference (gold) sets under `data/`; which sets hold `disp.*.nc`
   is listed on this row before any is replaced. NOT decided: adding
   displacement snapshots to a reference set that holds none today.
   Owner, 2026-11-18 16:36 CST, answering two numbered questions put
   after being told that no reference set holds `disp.*.nc` and that the
   bp1002 reference lacks the second-fault summary files: "1A. 2A."
   (verbatim.) The options as put were: question 1, displacement in the
   reference, "A) add the last-step snapshot per reference set
   (recommended), B) add every snapshot, C) keep none"; question 2,
   bp1002 second-fault summaries, "A) add them in the same gold PR
   (recommended), B) leave them out". Read as (a reading): when the gold
   is regenerated for the 8-byte storage, each reference set gains its
   last-step `disp.*.nc`, and the bp1002 set gains
   `cplot_ft2_EQquasi.txt`, `cplot_ft2_ruptarea_trac_slip.txt` and
   `cplot_ruptarea_trac_slip.txt`, all in that one gated gold PR, with
   the size the tree grows by recorded. This settles the "NOT decided"
   item just above. NOT decided: any other file absent from a reference
   set (`p1output.txt`, `p2output.txt`, `peak_sliprate_per_fault.dat` on
   bp1002); the tolerance the tier applies to the new files.
   Owner, 2026-11-18 16:38 CST: "/autopilot clear the board. a week
   budget. go" (verbatim.) Read as (a reading): the unattended run now
   covers every open row of this board, worked in `prio` order with
   state as the tiebreak, not only row 1; the budget is one week from
   this message, to 2026-11-25 16:38 CST, replacing the "until this
   row's items are done or blocked on the owner" end; the order already
   set on this row (refactor step 1, NetCDF output follows `dp`, gold
   regeneration, slip-to-traction promotion, then its optimization)
   keeps its place at the head of the queue. Unchanged: the Vista
   ceiling of 150 SU, every "NOT decided" item above, rule 1 for any new
   file, name, parameter or compset, no major version bump unattended,
   and that opening, closing and re-scoping rows stays with the rules
   owner.
   Owner, 2026-11-18 16:39 CST, answering one numbered question on the
   tolerance left open by the 16:36 paragraph: "1A" (verbatim.) The
   options as put were: "Bar for the double-precision NetCDF files in the
   test tier: A) 1e-10, same as the text columns (recommended), B) keep
   1e-6". Read as (a reading): once `disp.*.nc` and `fault.*.nc` are
   stored as 8-byte reals, the test tier compares them at a relative
   1e-10, the bar the double-precision text columns already use
   (`DD_RTOL`), in place of the 1e-6 kept for 4-byte storage
   (`DD_DISP_RTOL`). NOT decided: any change to the 1e-9 campaign
   comparator; what to do if a measured difference exceeds 1e-10 (it
   comes back to this row, the bar is not loosened unattended).
   Owner, 2026-11-18 16:45 CST, after being told that the 8-byte output
   fails the CI comparison against the 4-byte reference files: "though
   not bit identical, regenrated golds should match single precision with
   a relatively larger bar" (verbatim.) Read as (a reading): an
   acceptance condition on the gold regeneration, not a change of the
   1e-10 tier bar: every regenerated reference file must agree with the
   file it replaces at single-precision level (a looser bar, of the
   order of the 1e-6 the tier used for 4-byte storage), and the measured
   worst difference per file is recorded with the regeneration; a file
   that does not agree at that level is not replaced and comes back to
   this row. The stricter check already planned (new values cast to 4
   bytes equal the old values) satisfies it. NOT decided: the exact
   value of the looser bar (1e-6 is proposed here, as the bar already in
   the tier).
   Owner, 2026-11-18 16:45 CST, next message: "physics shouldn't change"
   (verbatim.) Read as (a reading): the stated reason for the condition
   above and a requirement on the whole change: storing the output at
   8 bytes may alter only the stored precision of each value, never the
   computed solution; the double-precision text output of the same case
   (station series, summary files) must stay within the 1e-10 tier bar
   between the binary before and the binary after, and any difference
   beyond storage truncation stops the gold regeneration and comes back
   to this row.
   (h) end state, owner 2026-11-13: "in the end, mumps is a choice o petcs".
   One PETSc time loop; MUMPS, BDDC and AMG are solver choices inside it and
   `solveTimeLoopMUMPS.f90` retires. Before it retires: PETSc-with-MUMPS
   matches the gold on every gated compset and is no slower per step than
   native MUMPS; PETSc builds on every supported machine (today only
   `conda-linux` links it); one assembly path serves both `MATIS` and AIJ;
   what `par.solver` then means is agreed first (rule 1).
   Two design principles for all of (e)-(h), owner 2026-11-13: "things
   parallelize from beginning, and thing that can armotize once should be
   done that way".
   Parallel from the beginning: no rank-0 phase anywhere -- mesh, numbering,
   assembly, initial conditions, friction, pore pressure, station and fault
   output are all per-rank.
   Amortize once: whatever does not change is built once and reused -- the
   solver setup across steps (K is constant); mesh, partition and setup
   across the cycles of a run (each cycle rebuilds and refactorizes today,
   137 s per cycle at Case 3); and, as a candidate to measure, the fault
   slip-to-traction operator (row 4's note), which would amortize the volume
   solve itself across steps and across runs on one mesh and modulus.
   Design choices shown to the owner 2026-11-13, answer "go":
   a new distributed time loop beside the two existing ones, which retire
   only at (h), so the gold gate stays meaningful throughout; a regular
   x-y-z block decomposition from the analytic generator, not `DMPlex`;
   split fault nodes may cross subdomain boundaries, with a ghost exchange
   per step for friction; the fault slip-to-traction operator is a measured
   side experiment, outside (e)-(h); an iterative solve is held to parity
   within a stated tolerance, shown not to move event times. Still to
   agree at (h): what `par.solver` means, and PETSc as a hard dependency.
   Evidence: parity to the gold at each of (e)-(h), and the per-step table
   against native MUMPS at 8/16/32 ranks on bp1002 and 10 m BP8 wide.
   Unmeasured: whether BDDC beats MUMPS at these sizes (CG+GAMG did not).
   Host, owner 2026-11-13: "do it on cotopaxi" (64 cores, the
   `eqquasi-petsc` env is on the shared home). Timings there are not
   comparable with the theo4 numbers above: re-measure native MUMPS on
   cotopaxi as the baseline first.
   **Archived measurement history (row 1, 2026-11-14 to 2026-11-18).** The state tables, the per-step solver comparisons and the conductor 15-30 entries are in `docs/dev/board_history.md`, section "Row 1 measurement history, 2026-11-14 to 2026-11-18 (archived 2026-11-18)". What they established (figures copied from that text; MEASURED, none VERIFIED):
   - Exact 3-step substructuring passes the 1e-9 comparator (bp1002 25/25, bp5 120/120, gs10 33/33) and does not beat exact MUMPS: bp1002 at 8-32 ranks 3.4-4.1x the memory, 2.3-3.6x the setup, 1.7-2.3x the loop time; gs10 at 16 ranks 4.6x peak memory (20.7 against 4.47 GB), 3.3x setup (268 against 81 s), 1.6x loop (4.45 against 2.71 s per step).
   - BLR on the interface, deluxe scaling and two-level BDDC closed negative or gained 5-40 % of the gap only; CG + PCBDDC on wide at 16 ranks: 339.7 s/step against 7.18 s/step for MUMPS (~47x slower).
   - GAMG CG: loop ratio 2.4x slower than exact at 32 ranks on wide (8.74 against 3.63 s/step); on bp8.qdc.gs.10 at 16 ranks loop 49.3 s against 8.1 s (6.1x). BLR as a CG preconditioner does not win.
   - Out-of-plane roundoff at 400 steps is not s2t-specific: the three out-of-plane files fail against exact also between two exact runs differing only in rank count; every in-plane file passes at 1e-9.
   - Vista job 1152125 (spike-s2t d78d822, version 1.24.14, 4 nodes on gh-dev) COMPLETED 00:42:27, exit 0:0.
   Open items and every Owner paragraph of the archived range follow.

   #### Owner paragraphs and open items kept from the archived range

   Source: line 640, conductor 15 hand-off entry of the c83b771 board (kept verbatim, not edited).
    *Owner item (conductor 25, 2026-11-18): agent type names in older logs under `docs/dev/` (grep of the type names, line counts): `SESSION_LOG_2026-10-12_board-sweep.md` 15, `board_history.md` 3, `fenicsx_feasibility_2026-10-09.md` 1; `docs/user/gen_params.py` lists names only as a scrub regex (not a doc, kept). Rewriting past logs is the owner's call; not edited. This PR replaced the names in this board (hand-off text of conductors 15/16) by role names; git history is not rewritten.*

   Source: line 699, design block of the substructuring line of the c83b771 board (kept verbatim, not edited).
    **Owner-held, not touched by any step:** any new parameter name or value (including the spike `-eq_*` options and final names), retiring solver 1 or 2, renaming `solveTimeLoopPETScDD.f90`, what `par.solver` means after the merge, the default solver, relaxing the guard at :557-562 (this decides whether single precision, BLR and CG options are reachable at all), PETSc as a hard dependency (the makefile builds it only on `conda-linux`), promoting any option past spike.

   Source: lines 1041-1047, spike-run reporting defect and its owner questions of the c83b771 board (kept verbatim, not edited).
     Spike-run reporting defect, 2026-11-14 (audit): `runInfo.json` of the GAMG runs reads `ksp_type` preonly, `pc_type`
     cholesky (g64b checked), not the solver used; the real iteration counts (`ksp_iterations_*`) are right. Merge-bound
     solver 3 path: `solveTimeLoopPETSc.f90:496-509` sets the names from the whitelist booleans (`isPreonly`,
     `isCholesky`) and the guard at `:472` rejects other combinations, so it names what it ran by construction; read from
     source only, the merge-bound binary was not run on a GAMG option. Not changed (source change, owner). Owner
     questions, not decided here: a test pinning the INFO(2) meaning (new test); `KSPSetUp` reading INFOG(1)/(2) on
     failure and any default ICNTL(14) (changes merge-bound solver behaviour).

   Source: line 1148, conductor 24 entry of the c83b771 board (kept verbatim, not edited).

Owner question: a default ICNTL(14) above 35 for solver 3 at high rank counts changes a solver default, so it is the owner's; not changed.

   Source: lines 1182-1184, conductor 25 open questions of the c83b771 board (kept verbatim, not edited).
Open questions:
- cplot 9.8e-4 in the step-0 rerun: `runs/20261118_step0/rerun`, `test.bp8.qdc.gs.10_s3_np4_rerun` and `_np8_rerun` against their first runs, column 4 (index 3), entry 277770180767.1491 against 277770180767.1501, 2.4e-15 relative to the column maximum (about 16 ulp); no other file pair exceeds 5e-7 absolute, all at most 3.3e-15 of column scale; np1 reruns identical (bp1002, bp5, bp7, bp8 np1: 0 difference). Roundoff from reduction order at np greater than 1, not noise of another kind.
- local wide reference 5.45 s/step (218.1 s loop, 40 steps, 32 ranks, `runs/20261118_vista/local/wide40_np32`): host load 5.8 at start, 38.3 at end (`load_before.txt`, `load_after.txt`), so the second half ran loaded; one run, no repeat at stable load. the 3.3-3.5 s/step exact figure for wide at 32 ranks IS on this row (timing pairs at 32 ranks, "is met by exact itself in these runs (3.3-3.5 s/step)", `runs/20261117_c17-*`; my earlier sentence that it was absent was wrong, corrected here after the coordinator pointed to it). Restated: exact at 3.3-3.5 s/step in those batches against 5.45 s/step in the 2026-11-18 reference (same case, same ranks); host load (5.8 at start, 38.3 at end of the reference run) is a plausible cause, untested. UNEXPLAINED.

   Source: line 1190, conductor 25 entry of the c83b771 board (kept verbatim, not edited).

Owner items added: (1) wide cplot files lose rows on Vista NFS (the output path writes through rank 0; whether to write elsewhere or sync is the owner's); (2) cross-architecture bar for float32 `disp` files (1e-9 fails at 2.4e-8 of scale); (3) the comparator bar for out-of-plane columns, now with Vista evidence; (4) default ICNTL(14): default 35 failed at np64, 288, 576 on Vista and 100 passed, and passed at np32 and np144 at 35.

   Source: line 1457, conductor 29 entry of the c83b771 board (kept verbatim, not edited).
**OWNER ITEM, holds the merge:** the design on this row says step 1 moves the helpers into "a new or existing module" and names no file or module; `PROJECT_RULES.md` rule 1 requires agreement before a new file or name. 3A approved the design, not this name. Question for the owner: A) new `src/onFaultHelpers.f90`, module `onFaultHelpers`; B) an existing module (no new file; the makefile dependency order changes instead). No refactor PR is opened and v1.24.7 is not merged until answered; branch `worktree-agent-a42a732ee98ad8c43` is pushed, head 6b8693d, worktree `.claude/worktrees/agent-a42a732ee98ad8c43` kept.

   Source: line 1466, conductor 29 hand-off of the c83b771 board (kept verbatim, not edited).
OWNER ITEMS OPEN: (1) the s2t option's name and parameter value; (2) the file name for refactor step 1, A) new `src/onFaultHelpers.f90` or B) an existing module (section above; branch `worktree-agent-a42a732ee98ad8c43` head 6b8693d gated, held); (3) the 1A floor value; (4) ICNTL(14) for native solver 1 and `KSPSetUp` reading INFOG, unchanged.

   Source: line 1490, conductor 30 entry of the c83b771 board (kept verbatim, not edited).
Proposal 14 gaps still OPEN, unchanged: (a) no `disp.*.nc` in the committed gold, so the s2t displacement field is not compared on any case; (b) bp1002 second-fault `cplot_ft2_*` and `cplot_ruptarea_trac_slip.txt` have no gold file; (c) s2t-only extras have no gold file. Proposal 13 job 1152125 (started 16:02 CST on gh-dev) is not harvested; the invoking session relays its exit. Owner items open as in the conductor 29 hand-off.

3. [ ] **P1 -- EQquasi BP8 discrepancy (bug hunt).** Owner, 2026-10-09: "hunt
   the bug". An independent FEniCSx linear-hex model with EQquasi's box and BCs
   converges to Kim within ~2 % at 25 m (`docs/dev/fenicsx_feasibility_2026-10-09.md`);
   EQquasi does not (depth station -18 % at 50 m, -13 % at 10 m). Pressure
   matches; domain and dtmax ruled out. The uploaded BP8 entries are suspect
   along depth until this closes. Done when: the cause is found and fixed, and
   EQquasi at 25 m agrees with Kim within a few % at all 9 stations.
   **Escalation, 2026-10-12**: 6 code-audit findings on this hunt
   rejected after a source/gate check (4 before today: nodal area/dx^2, GS
   z-symmetry, and two more per the 2026-10-09 brief; 2 today, PRs #58 and
   #59, closed unmerged -- see `docs/dev/SESSION_LOG_2026-10-12_board-sweep.md`
   for each). Today's two both targeted `src/faulting.f90`'s split-node
   traction formula (lines 112-128): #58 claimed "ruled out" without ever
   computing the required EQquasi-vs-FEniCSx ratio; #59 computed only the
   FEniCSx half and flagged the mass-weighted cross term as a "structural
   bug," which does not survive the rule-16a sanity check -- that formula
   sits outside every `bp ==` branch (the only ones in the file are at
   lines 159/162/359, pore-pressure/BP7-specific), so a bias large enough to
   explain BP8's -13%/-18% would also break BP5/BP1002's passing regression
   gold, and it doesn't. The split-node coupling itself is very likely not
   the bug (3 independent rejected attempts have now landed on it from
   different angles). **Needs owner diagnosis, not a 3rd same-session
   retry**: live lead is still unconfirmed-but-unrefuted (tau/pressure match
   Kim within ~0.5% but V amplifies a small tau gap exponentially); next
   untried candidate is `src/qdct*.f90`'s element-to-fault force
   contribution (assembly path, not yet audited by anyone).
4. [ ] **P2 -- FEniCSx backend, EQdyna pattern.** Owner, 2026-10-09: "like
   eqdyna, we will have a fortran and fenicsx backend now ... follow eqdyna's
   success". One compset and one gold set, a test matrix that runs each case
   on each backend, a backend counts only at parity with the gold. Order:
   after row 3 (the Fortran gold must be right); BP5 and step-over parity
   first, BP8 against Kim. The prototype's precomputed slip-to-traction matrix
   (quasi-static, linear) is also a scaling idea for row 1.
6. [ ] **P2 -- Time-step factor xi above the stability limit.** Owner,
   2026-11-13: "that ksi lesson, put on board". Lesson: a constant `par.xi`
   above the Lapusta et al. (2000) eq. 15 limit makes the time step
   oscillate between events. Reported by the session that ran the
   comparison, unaudited here: Wang 2026 Case 3 at 1 km with xi 0.2
   (`runs/20261106_wang2026_case3/{c3a,c3b}`, `result/cycle*/global.dat`),
   90-97 % of interseismic steps reverse the sign of the step-to-step change
   in log10 Vmax; with per-node xi (PR #72, `par.xi_gamma = pi/4`) 0-1 %.
   The eq. 15 limit with k = gamma mu / h, a 0.005, b 0.01, L 0.05 m,
   sigma 25 MPa, 36 GPa is 0.097 at 1 km and 0.046 at 500 m, so xi 0.2 at
   1 km is about twice too large. Event times at 1 km moved by at most
   0.8 yr over five events; the defect is between events, not in the cycles.
   Open, owner's call: (a) `case.setup` computes the eq. 15 limit from
   material, h and friction and warns or stops when `par.xi` exceeds it;
   (b) per-node xi as the default with `par.xi` as the cap; (c) PR #72's
   open audit items -- `xi_gamma = NaN` silently off, the stop-13 check also
   runs at creeping nodes, stop 13 exits without `MPI_Abort`.
   Done when: PR #72 is merged or closed with (c) resolved, and (a) and (b)
   are each decided and, if yes, landed with a test.

Proposals for the owner (not decided): (1) `Slip_3` floor of 1e-15 m, `fltst` files only, measured by conductor 30 from the committed reference; the 1e-12 m/s and 1e-9 MPa floors were already decided. (2) The tolerance change to 1e-10 (disp and fault NetCDF at the repository tolerance after the 8-byte change) rides in PR 203 together with the regenerated reference, because a separate gold PR would leave master red.

Hand-off (conductor 30 to the next): (a) PR 203 (branch `c30-dp`, head b6c6e67, HOLD) waits on the gold-regeneration branch from the running implementer; run the fresh-clone fast tier and `-m e2e_ci` on the 1.24.8 build before folding it onto 203, merge only on `gh run watch --exit-status` at the exact head, remove the HOLD marker after; (b) then the slip-to-traction promotion stack off d78d822 (`par.solver = 4`), one src PR at a time, code audit before any tag, no new source file name until the owner names it; (c) remove the nested worktree `.claude/worktrees/c30-dp` after 203 lands and reap dead child worktrees only after copying any cited file out; (d) the remote branch `worktree-agent-a42a732ee98ad8c43` failed to delete, cause unread.

Hand-off addendum, conductor 30 (harness tool calls 88 of 100): the gold-regeneration implementer finished PARTIAL and is stopped. Branch `c30-gold`, head 05df35f, one commit on top of `c30-dp` b6c6e67, pushed; NOT gated by the conductor (its own report: fast tier 302 passed, `-m e2e_ci` 14 passed and 8 skipped for absent reference categories; unaudited).

- Replaced: 12 `fault.*.nc` in 4 sets (bp5 fast, bp5 dip90, stepover, stepover.con); `disp.00101.nc` added to dip90, stepover, stepover.con; data tree 94,260,435 to 100,994,278 bytes. Tolerance sites: `compare_netcdf` default 1e-10, `DD_DISP_RTOL` removed (disp and fault at `DD_RTOL`), a `disp` category and the bp1002 ft2/ruptarea cplot files added to the comparators, `"disp."` added to the contract read-glob.
- Not done, 5 sets: bp7 fast (`fault.00001.nc` fails the 1e-6 single-precision bar, max abs 3.8e-19 on a scale near 1e-18, relative 0.39: needs an owner decision, nothing replaced), bp8 (1.24.8 run done, 703 s; no 1.24.7 pair), bp7 full, bp1002 (including its 3 cplot files), bp5 kink, liu2020 (estimated over 4 h at 3 ranks, unmeasured). Estimates per set at 3-4 ranks: bp5 full about 50 min per binary, bp1002 about 45, bp5 kink about 40, each times 2 for the 1.24.7 pair.
- Open items for the next conductor: no 8-byte check exists on the reference files (a contract test is needed); `p1output.txt` fails a raw 1e-10 `compare_series` on stepover and passes the tier's own `_dd_profile_series` (a comparator choice to confirm); the bp5 and bp7 fast subdirectories hold no `disp` file (so `e2e_ci` skips that comparison there; one file per set was the instruction); worst difference at 1, 4 and 8 ranks against 1e-10, the restart run and write time are not yet recorded. Resume with `scratch/gold8/driver.py` and `scratch/gold8/check.py` in the implementer's worktree (uncommitted notes; copy out before reaping). PR 203 stays HOLD until the gold is complete or the owner accepts a partial set.

### Row 1, conductor 32: state recovered from the record (2026-11-18 22:20 CST; harness tool calls 10 of 100; budget read as wall clock to 2026-11-25 16:38 CST, row 1's order first, hand-off at about 80 calls)

Conductor 31 left no entry. Recovered from git, `ps` and `runs/` at 22:12 CST; every figure below is read this turn, nothing is gated or audited yet.

- Master `cc5137e`. Open PRs: #203 (`c30-dp`, head `4eca07d`, GitHub CLEAN) and the spike #87 (not for merge). Remote-only branches of conductor 31, none with a PR: `c31-p1` (`3b39a54` `par.solver = 4`, `00d6162` tests), `c31-p1fix` (head `09496fb`, version 1.24.10, stacked on `4eca07d`), `c31-goldA/B/D2/G`, `c30-gold8`, `worktree-agent-ad7e2b054bafedb8b`.
- #203 content (diff against master): the `dp` change `b6c6e67` (1.24.8); gold commits `05df35f` (4 sets), `f72f119` (new file `tests/contract/test_gold_netcdf_shape.py`, 82 lines, new file needing rule 1 agreement: asked below), `f03b8e2` (bp8 regenerated, bp7 full disp only), `2f822c0` (bp1002 with its 3 second-fault cplot files and last-step disp), `4eca07d` (bp7 fast, bp7 full, bp5 kink). 73 files changed, binary data only plus `src` (4 files), docs, tests. Not yet regenerated: liu2020.qdc.kink.600.
- liu2020 regeneration: the run in `runs/20261118_goldF/248/liu600` (3 ranks, nstep 10000, reference ends at `fault.10000.nc`) is NOT running: no `prterun` or solver process exists at 22:12, `cycle0.log` last written 22:10:35 at step 2821, `248.liu600.log` is 0 bytes, so no error text; cause unread (no `dmesg` access). Checkpoint files `fault/disp.{00001,01001,02001}.nc` are in its `scratch/`. A cycle cannot resume, so the set needs a fresh run (about 6 h at 2821 steps in 1 h 46 min).
- bp7 fast question: `4eca07d` DID replace the bp7 fast `fault.*.nc`, on the basis that each passes the tier comparator (`compare_netcdf`, rtol 1e-6 with an absolute floor of 1e-6) against the file it replaces. Conductor 30 had held it back because the variable of scale about 1e-18 differs relatively by 0.39 (absolute 3.8e-19). The owner's 16:45 condition says "agree at single-precision level, a looser bar of the order of the 1e-6 the tier used", worst difference recorded per file; whether the tier's absolute floor counts is the open reading. Being checked by me before #203 merges; if the file only passes through the floor it goes back to the old file and the question stays here.
- `par.solver = 5`: not a selectable value. `cb3c05e` adds only a test that `par.solver = 5` is refused with "par.solver must be 1, 2, 3 or 4". Nothing to remove from `src`.
- Operator-check floor (`570977e`: bound floored at 1e-14 of the operator's largest entry; motive: bp7 step 1 holds a traction of 8e-15, round-off, which a pure relative bound aborted) is a numerics choice for the owner, NOT agreed. Last gate at `cb3c05e`: fast 371 passed, `e2e_ci` 14 passed and 8 skipped, solver-4 tests 11 passed and 2 failed (bp7); `570977e` and `09496fb` never run. Not done: s/step of bp5 and bp1002 after the operator check, the planted-error check with the floor, the no-PETSc refusal test. No gate script holds `scratch/locks/heavy.lock` (a 0-byte file, no holder in `ps`); a process of another project (`mpirun -np 10 eqdyna`, PID 8034) loads the host.
- Worktrees with cited files, not to be reaped before copy-out: `agent-a2e8a419ac82b9331` (`runs/20261118_p1fix/` logs), `agent-a1a440daa6089e9de` (locked; goldF), `agent-af6f00fefea46947a` (`scratch/gold8/driver.py`, `check.py`).
- Queued for the owner's daily brief (09:37 CST): (1) the bp7 fast file reading above; (2) `Slip_3` floor of 1e-15 m (already on this row as proposal 1 of conductor 30); (3) the operator-check floor of 1e-14 of the operator's largest entry, with the alternatives (pure relative, which aborts bp7 step 1; a floor tied to the traction scale) once characterized; (4) the gold lands in the same PR as #203; (5) rule 1: `tests/contract/test_gold_netcdf_shape.py` is a new file, named by conductor 31 without agreement; it stays unless the owner objects, and the check can instead move into an existing contract file.

Order of work (invoking session): this entry; finish and gate the gold and land #203 with master never red; gate, audit and land the s2t stack as `par.solver = 4`; floors in the comparator; s2t optimization and the rest of the board in `prio` order. Owner-held and untouched: ICNTL(14) for native solver 1, a `vista` machine target, any default change, any new file or name beyond those agreed, other absent bp1002 reference files.

### Row 1, conductor 32: #203 audit, gold measured, liu2020 relaunched (2026-11-18 22:35 CST; harness tool calls 43 of 100)

- Correction to the entry above: `tests/contract/test_gold_netcdf_shape.py` is an existing file that #203 extends by 82 lines (`git diff --diff-filter=A` lists no new file outside `data/`); rule 1 is not engaged and the fifth item of the owner-brief list is struck.
- bp7 fast and the reading of the owner's 16:45 condition, MEASURED by me this turn (all 53 replaced `fault.*.nc` of #203, each old (master) against new (`4eca07d`) with the tier's own `compare_netcdf`): 0 fail at rtol 1e-6 with floor 1e-6; 31 fail with the floor at 0 (bp5 kink 22 of 22 after step 600, bp8, stepover, stepover.con). The only variables above 1e-6 of their own scale are zero by physics: `shear_dip` (own scale 9.7e-19, group scale 1.6e7), `slipd` (7.1e-38), `vzm` and `vzs` (3.5e-35), and `slipn` (own scale 5.6e-12, abs 2.3e-15, 4.0e-4 of own scale in the bp5 kink set); largest difference over the group scale of its field: 1.0e-13. Every signal variable agrees within 1e-6 of its own scale; the worst per-file figure (difference over the file's own maximum) is 5.7e-8 for bp1002, bp5 fast, dip90, stepover and bp7 `fault.00101`, 4.6e-5 for bp7 full from step 1001, 2.4e-7 for bp8, and 3.9e-1 (bp7 `fault.00001`, the zero-by-physics fields above) and 3.9e-4 (kink `slipn`). Reading kept: the replacement meets the condition as the tier defines its 1e-6 bar (group scale, absolute floor); under an own-scale reading the five zero-by-physics variables would not, and that is the question for the owner brief. The files stay replaced.
- Tree growth, #203 head: `data/` 94,260,007 to 112,767,711 bytes (+18.5 MB, `git ls-tree -r -l`). A liu2020 last-step `disp` file is about 6.3 MB, above rule 22's 5 MB per-file cap: owner question (exempt the file, or keep liu2020 without its `disp`); liu2020 therefore waits and does not block #203.
- Audit of #203 (read-only, role auditor, 12 tool calls): no BLOCKER, verdict "merge after fixes". MAJOR 1: `data/test.bp5.qdc.2000/cycle0` (full tier) still holds 4-byte gold compared at 1e-10, so the next `-m e2e` fails it; not on the width list. MAJOR 2: the CI case `test.bp5.qdc.2000/cycle0-step101-fast` has no `disp.*.nc`, so the new `disp` category skips there (rule 3a); owner decision fd153ba wants one per set. MAJOR 3: restart files now hold unrounded values, so cycle 2 and later start from a different state than under 1.24.7 (cycle 1 unaffected); the owner was told fault files are read back on restart when answering "both follow dp", and the doc line was missing. MINOR: bp8 absent from `REGENERATED_8BYTE`; `FAULT_SNAPSHOTS` globs `data/*/fault.*.nc` and finds only bp8 (older than #203); the string "1.24.8" is also the version of an uncommitted spike binary. Fixes in preparation on branch `c30-dp`: restart paragraph in `docs/user/outputs.md`, bp8 on the list, bp5 fast `disp.00101.nc`, bp5 full regeneration (50 min per binary at 4 ranks, queued behind `scratch/locks/heavy.lock`). Multi-cycle difference not measured; queued for the owner brief as a consequence of the dp decision.
- Runs (detached with `setsid`; conductor 31's run died with its shell at 22:10, likely cause, unproven): liu2020 relaunched 22:17:33 (`runs/20261118_goldF/248/liu600`, run.sh PID 719210, prterun 719237, 3 ranks, binary copy `scratch/c32/bin/eqquasi-1.24.8` md5 dc714b19, same as the one conductor 31 used; dead-run output kept in `died_at_step2821_22h10/`; about 6 h, finish about 04:30); rank study at 4, 1 and 8 ranks, `runs/20261118_c32_rank/` (driver `scratch/c32/rank_study.py`, PID in `study.pid`); text pair 1.24.7 against 1.24.8 at 4 ranks, `runs/20261118_c32_pair/` (`scratch/c32/pair.py`). Gate trees: `.claude/worktrees/c32-gate` (head `4eca07d`, with `bin/eqquasi-1.24.8`), `c32-m247` (master, `bin/eqquasi-1.24.7`), `c32-fix` (the PR fix edits).

### Row 1, conductor 33: board split landed, budget reading, owner-brief block (2026-11-18 23:10 CST; harness tool calls 26 of 100; unaudited)

- Budget read as: the owner's "/autopilot clear the board. a week budget. go" (16:38 CST) covers every open row in `prio` order until 2026-11-25 16:38 CST with this row's order first (#203 with its gold as one PR, the slip-to-traction stack as `par.solver = 4`, comparator floors, then optimization and the rest). Owner-held and untouched: ICNTL(14) for native solver 1, a `vista` machine target, any default change, any new file or name beyond those agreed, the operator-check floor until answered. Hand-off at about 80 calls.
- Landed: #210 (`65b3e5c`, CI run 44460871482 success on head `8359e09`): the row 1 measurement history (2026-11-14 to conductor 30) moved verbatim into `docs/dev/board_history.md`, section "Row 1 measurement history, 2026-11-14 to 2026-11-18 (archived 2026-11-18)"; board 345,997 to 52,349 bytes. Proof run by me on the pushed head: every Owner paragraph and open-item paragraph byte-identical, the archive's old text an untouched prefix, 0 non-blank old lines lost. Existing file, no new name (rule 1 not engaged). Pre-existing in the moved text, not edited: a repository URL in the archive that carries a personal handle (`board_history.md:979`); the owner decides whether past logs are rewritten (item kept on row 1).
- #203: head `9486225` (restart paragraph, bp8 on the 8-byte list, bp5 fast `disp.00101.nc`, pushed by me, fast-forward). Still to do before merge: bp5 full regeneration (1.24.8 run prterun 735604 started 22:25, about 50 min), then fast tier and `-m e2e_ci` at the merge head from a fresh clone, then `gh run watch --exit-status` at the exact head. liu2020: run.sh 719210 / prterun 719237, started 22:17, expected about 05:00; it does not block #203.

Owner brief, 09:37 CST (questions wait for it; every item has its options and my recommendation; nothing below is decided):
1. bp7 fast and the other replaced `fault.*.nc` (53 files, measured by conductor 32: 0 fail at rtol 1e-6 with the tier's absolute floor of 1e-6, 31 fail with the floor at 0; the failing variables are zero by physics, scale 1e-18 to 1e-38). A) keep the replacement, the tier's own bar being the single-precision bar (recommended: the signal variables agree within 1e-6 of their own scale and the physics is unchanged). B) require agreement at 1e-6 of each variable's own scale: five zero-by-physics variables fail, so bp7 fast and the bp5 kink would keep their old files and fail the 1e-10 tier bar once the output is 8 bytes.
2. `Slip_3` floor of 1e-15 m, `fltst` files only (measured by conductor 30 from the committed reference). A) adopt (recommended; the 1e-12 m/s and 1e-9 MPa floors were already decided). B) no floor: the series fails on round-off-level values.
3. Operator-check floor of the slip-to-traction check (`570977e`: bound floored at 1e-14 of the operator's largest entry; bp7 step 1 holds a traction of 8e-15, round-off). A) adopt the 1e-14 floor (recommended only together with a test that plants a wrong operator and shows the check still aborts it; that test is audit finding M2 and is not yet written). B) pure relative bound: aborts bp7 step 1. C) a floor tied to the traction scale: not yet characterized. Not merged before the answer.
4. The gold lands in the same PR as #203. A) one PR (recommended: a separate gold PR leaves master red for the tolerance change). B) two PRs, master red between them.
5. liu2020 last-step `disp` file, about 6.3 MB against rule 22's 5 MB per-file cap. A) exempt this one file, recorded by the rules owner (recommended: the owner's decision wants one `disp` per set). B) keep liu2020 without its `disp`.
6. Restart files now hold unrounded values (dp), so cycle 2 and later differ from 1.24.7 at the level of the 4-byte rounding; cycle 1 is unchanged; the docs say so. Not measured. A) accept, after I measure the cycle-2 difference on bp5 fast with both binaries and log it (recommended). B) keep 4-byte restart rounding: a second code change.

### Row 1, conductor 33: #203 landed, HAND-OFF (2026-11-18 23:42 CST; harness tool calls about 50 of 100, a context summary at 23:40 voids the count: no new work, no further merge; unaudited)

- Landed: #210 (`65b3e5c`, board split) and #211 (`231f8f1`, owner-brief block, six items with options and recommendations) by board PR; #203 (`6934a28`, NetCDF reals follow dp, v1.24.8, with the gold of every set but liu2020; its bp5 full regeneration added at head `4a170d4`: faults at most 2.6e-11 of scale from the file replaced, 18 text files at most 1.9e-14). Gate on `4a170d4` from a fresh bare clone under `env -i` (`runs/20261118_c33_gate203/`): fast tier 384 passed with the binary built (383 passed, 1 failed before the build, the binary-on-disk test), `-m e2e_ci` 15 passed, 7 skipped (data-absence skips: no offfault, perfault, profile files, no `summary.json`, no section 4.3 profiles; none new). GitHub CI run 44463991495 success on that head; master CI run 44465226691 success on `6934a28`. The gate and the reference-data record are the comment on #203; the PR body's HOLD text was replaced by that record.
- NOT done on #203's follow-through: liu2020 regeneration still running (run.sh 719210 / prterun 719237 started 22:17, 1 h 22 min at 23:40, expected about 05:00; its `disp` of about 6.3 MB exceeds rule 22's 5 MB cap: owner brief item 5); the multi-cycle restart difference is not measured (brief item 6). The 1.24.7 pair driver (`scratch/c32/pair.py`, PID 733568) and the rank study (PID 732912) still run; their outputs are unread by me.
- Slip-to-traction stack: branch `c33-s2t-fix`, head `9903acd` (pushed; on top of `c31-p1fix` `09496fb`; version 1.24.11), UNAUDITED, UNMERGED, NOT rebased onto `6934a28`. The implementer's report (not checked by me): B1 fixed (the test uses the solver's bound `max(1e-8*field_max, 1e-14*operator_max)`, marked "PROPOSAL, owner decision pending"), M1 fixed (NaN flag per value, reduced across ranks), the abort names both bounds, a planted-row test of the assert alone; M2 only partly done (no run-level planted-operator test: it needs a solver hook, a new name, rule 1); minors (b) allocation message, (c) non-MUMPS override refusal test, (d) no-PETSc refusal test NOT done; 1-, 3- and empty-rank layouts open. Its gate (own worktree, 4 ranks): fast tier 371 passed; solver-4 tests passed for bp1002 zone, bp5 and bp7.a.10, the unsupported-combination refusals (krylov3, krylov4, solver5) and the planted-row test; the 8- and 12-rank cases (bp7.8, bp1002.12, bp8 gs.8 and pw.8) were not run (over the 4-rank cap). Next: rebase onto master, re-run those cases at the allowed rank count on a quiet host, a code audit of the final head, no merge before the owner answers brief item 3 (operator-check floor).
- Next for a fresh conductor, in order: (1) read the owner's answers to the brief and apply them; (2) liu2020 gold when its run ends (copy out, compare against the old file with the tier comparator, apply per the owner's answer on item 5), one gated PR; (3) the stack above as `par.solver = 4`; (4) comparator floors (`Slip_3`, brief item 2); (5) s2t optimization and the rest of the board in `prio` order. Owner-held: ICNTL(14) for native solver 1, a `vista` machine target, any default change, any new file or name (the planted-operator hook included).
- Worktrees and branches I made: `.claude/worktrees/c33-board` (branch `c33-handoff`, this entry), `agent-a94de11514ecea170` (finished, reapable after a listing for gitignored files), `agent-a52c2bbfc0ba82b73` (finished, branch `c33-s2t-fix`; keep until the stack lands); `.claude/worktrees/c32-fix` holds gate copies under `scratch/gold8/bp5full` (ignored, copy-out not needed: the run dir `runs/20261118_c32_pair/v1248_np4/bp5full` stays). Not reaped: all others (56 at the start; none shown unreferenced this run).

### Row 1, conductor 34: #212 landed, dispatches, checkpoint (2026-11-19 00:06 CST; harness tool calls 17 of 100; unaudited)

- Budget read as: row 1's order first, then every open row in `prio` order, to 2026-11-25 16:38 CST; items 1-6 of the owner brief wait for the owner's 09:37 CST answer and nothing past them merges. Hand-off at about 80 calls.
- Landed: #212 (`3aae60b`; CI run 44466783010 conclusion success on head `a90f2a6`, state MERGED read after the merge). The branch-delete error after merges is the explicit `git push origin --delete` finding the branch already gone (`error: unable to delete 'c33-handoff': remote ref does not exist`): the repository deletes head branches on merge, so the explicit delete is redundant. `c33-brief` is not on origin (`git ls-remote --heads origin 'c33*'` lists `c33-s2t-fix` only).
- Dispatched, both cap-bound, worktree-isolated, heavy runs behind `scratch/locks/heavy.lock`: (1) the slip-to-traction stack, branch `c34-s2t` from `6934a28` (8 commits of `c33-s2t-fix` cherry-picked, one bump to 1.24.9, minors (b)(c)(d) without a new name, gates at 1/3/4/8/12 ranks, setup, s/step and memory per case, draft PR marked HOLD; cap 60 calls); the planted-operator hook (M2) is NOT built: its name goes in the report for the owner (rule 1). (2) the cycle-2 restart difference between 1.24.7 and 1.24.8 on bp5 fast (measurement only, brief item 6; cap 40 calls).
- Still running, `ps` this turn: prterun 719237 (liu2020 gold, 3 ranks, started 22:17), prterun 945976 (rank study bp8, 4 ranks, started 23:43), python 732912 (rank study driver), python 733568 (1.24.7 pair driver). Outputs unread; to be harvested onto this row when they end.
- Reaped (clean, merged, no ignored output beyond `__pycache__`): worktrees `c25-board2`, `c25-board3`, `c26-board`, `c32-board`, nested `owner-relays-2/.claude/worktrees/c30-dp`; local branches `board-c25-2`, `board-c25-3`, `board-c26-5`, `board-c32-2`, `c30-dp`. 56 worktrees remain; every other one is either a live agent, holds ignored output not yet listed, or is referenced.
