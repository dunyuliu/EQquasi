"""The e2e case table and the one runner that executes any of them.

There used to be a file per benchmark -- five of them, ~950 lines -- each with
its own `create.newcase -> case.setup -> run.sh -> compare` fixture differing
only in a compset name, a few parameter overrides, and which reference to diff
against. Two consequences, both of which bit:

  - a change to the workflow had to be chased through every file. The Q0 ->
    cycle0 rename broke `tests/e2e/test_benchmarks.py` because it was updated in one place and
    not another.
  - each file grew its own comparison logic, so BP8 got full-array diffs while
    BP5 got none, for no reason other than which file was edited last.

Every case produces the same families of file:

    cycle0/  fault.NNNNN.nc  disp.NNNNN.nc  fault.r.nc  disp.r.nc
             fltst_strk*     srfst_strk*    global.dat
             cplot_EQquasi.txt  tdyna.txt   runInfo.json

and every difference between benchmarks is readable from the data rather than
from the benchmark's name: station column count (9 for BP5/BP7, 11 for BP8),
extension (.txt against .dat), `nid_fault` size (1, or 2 for the step-over),
whether the section 4.3 profiles exist (BP8 only), and whether anything
ruptured. So the comparisons are driven by what the reference directory
contains, exactly as `scripts/plotStations.py` decides its layout.
"""

import bisect
import glob
import json
import math
import os
import re
import shutil
import subprocess
import sys
import time

import pytest

from conftest import ROOT

sys.path.insert(0, os.path.join(str(ROOT), "scripts"))

# name, compset, parameter overrides, reference subdirectory, tier.
#
# "fast" cases also carry the e2e_fast marker and are what CI runs on every
# push. "full" cases run only under -m e2e. Costs measured on a 64-core shared
# host, 4 MPI ranks: 101-step cases ~90 s, BP8's 30 days ~800 s, a full BP5
# cycle ~3100 s.
CASES = [
    # fast: what CI runs on every push.
    ("test.bp5.qdc.2000", "test.bp5.qdc.2000",       {}, "cycle0-step101-fast", "fast"),
    ("test.bp7.qdc.a.10", "test.bp7.qdc.a.10",       {}, "cycle0-step101-fast", "fast"),
    # HPC_ncpu is overridden to 2, matching bp5 and bp7. The compset asks for
    # 20, which is right for the machine it was written for and impossible on
    # a GitHub runner: ubuntu-latest has 4 cores, so `mpirun -np 20` fails
    # before the solver starts. BP8 failed in CI for several releases for this
    # reason alone, passing locally every time on a 64-core host, and the
    # failure was unreadable because run.sh's log was never surfaced.
    ("test.bp8.qdc.gs.10", "test.bp8.qdc.gs.10",
     {"xi": 0.2, "nstep": 8000, "nt_out": 8000, "HPC_ncpu": 2},
     "",                                                            "fast"),
    # De-orphaned 2026-08-15 (references frozen the same day, UNVERIFIED
    # numbers, regression locks). stepover is the only ntotft > 1 row cheap
    # enough for the fast tier: it closes the rule-12 gap by putting
    # per-fault on-fault-input routing into every-push CI.
    ("test.bp5.qdc.dip90.2000", "test.bp5.qdc.dip90.2000", {}, "cycle0",               "fast"),
    ("test.stepover.qdc.1000", "test.stepover.qdc.1000",  {}, "cycle0",               "fast"),
    # Constraining twin: one flip (far_vel_load sign). Same cost.
    ("test.stepover.qdc.con.1000", "test.stepover.qdc.con.1000", {}, "cycle0",          "fast"),

    # The big ntotft > 1 row (the stepover fast row above is the cheap one):
    # three multi-fault bugs in this project's history were found by hand
    # because nothing in the gate exercised per-fault routing of on-fault
    # input (PROJECT_RULES rule 7).
    #
    # BP1002 is deliberately NOT in the fast tier. Until 2026-08-15 it was the
    # only ntotft > 1 case anywhere, so per-fault routing of on-fault input
    # was absent from every-push CI (the rule-12 gap); the stepover fast row
    # above now covers it. BP1002 still exercises the bigger mesh.
    #
    # The reason is cost: 41472 elements, 3821 steps, ~2600 s on 3 ranks of a
    # 64-core host. A GitHub runner has four cores. It stays in the full tier
    # until either the mesh is cheaper or CI runs somewhere that can afford
    # it.

    # full: -m e2e only.
    ("test.bp5.qdc.2000", "test.bp5.qdc.2000",
     {"nstep": 100000, "nt_out": 1000},          "cycle0",         "full"),
    ("test.bp7.qdc.a.10", "test.bp7.qdc.a.10",
     {"nstep": 100000, "nt_out": 1000},          "cycle0",         "full"),
    ("bp1002.qdc.2500", "bp1002.qdc.2500",   {},       "cycle0",         "full"),
    # Row-7 zoning experiment (PATHWAY_FORWARD, owner decision 2026-09-24):
    # bp1002 with the VW patch held fixed in area and distance from the
    # step-over while the segments are made unequal (A 57.5 km, B 102.5 km).
    # UNVERIFIED: no reference is committed yet -- the case is a first read of
    # a genuinely new experiment, not a regression lock -- so _run() skips it
    # ("no reference at ...") until data/test.bp1002.qdc.zone.2500 is
    # frozen from a reviewed run. 155 km mesh, ~1.3x bp1002's cost.
    ("test.bp1002.qdc.zone.2500", "test.bp1002.qdc.zone.2500", {}, "cycle0", "full"),
    # BP5 dip90 with only the surface kinked 10 deg (user-designed control;
    # the result it locks: rupture crosses the bend freely under BP5
    # friction). ~40 min on 3 ranks -- full tier.
    ("bp5.qdc.kink.2000", "bp5.qdc.kink.2000", {},       "cycle0",         "full"),

    # Recreated later, through this mechanism, from a clean run: BP5-dip90.
    # Its old reference was produced ad hoc, before the workflow rule and
    # before the gold/ flatten, so it is not carried forward -- a reference
    # nobody can regenerate is worse than none. Add a row here when its
    # reference exists and has been verified.
]

RUN_TIMEOUT_S = 7200


def case_id(c):
    name, _, over, ref, tier = c
    return f"{name}.{ref}" if ref else name


def reference_dir(name, ref):
    d = os.path.join(str(ROOT), "data", name)
    return os.path.join(d, ref) if ref else d


def env():
    e = dict(os.environ)
    e["EQQUASIROOT"] = str(ROOT)
    # bin/ first: a stale EQquasi elsewhere on PATH shadows this one and dies
    # with a shared-library error that says nothing about the real cause.
    e["PATH"] = f"{ROOT}/bin:{ROOT}/scripts:" + e.get("PATH", "")
    e["OMP_NUM_THREADS"] = "1"
    return e


def apply_overrides(udp, over):
    """Rewrite `par.<key> = ...` lines in place, keeping everything else."""
    if not over:
        return
    out = []
    for ln in open(udp):
        key = ln.split("=", 1)[0].strip() if "=" in ln else ""
        if key.startswith("par.") and key[4:] in over:
            out.append(f"{key} = {over[key[4:]]}\n")
        else:
            out.append(ln)
    open(udp, "w").writelines(out)


def installed_binary():
    """Path to bin/eqquasi-<the version src declares>, or None.

    bin/ holds versioned binaries and deliberately no plain `eqquasi`: one
    unversioned slot cannot hold two builds, and a run that spans a rebuild
    picked up the new binary on its next cycle.
    """
    v = declared_version()
    if not v:
        return None
    p = os.path.join(str(ROOT), "bin", f"eqquasi-{v}")
    return p if os.path.exists(p) else None


def declared_version():
    """The version src/globalvar.f90 declares, or None if it says nothing."""
    src = os.path.join(str(ROOT), "src", "globalvar.f90")
    if not os.path.exists(src):
        return None
    for line in open(src):
        if "EQQUASI_VERSION" in line and "'" in line:
            return line.split("'")[1]
    return None


def binary_version(exe):
    """The version bin/eqquasi prints in its banner.

    Raises rather than returning None when the binary cannot be run or prints
    no version: a check whose failure mode is "silently pass" is worse than no
    check, and this one exists precisely because a stale binary is otherwise
    undetectable.
    """
    r = subprocess.run([exe], cwd=str(ROOT), capture_output=True,
                       text=True, timeout=120)
    out = r.stdout + r.stderr
    # Anchored on the banner text, not on "any dotted number": a build date
    # like 2026.08.12 in the output otherwise reads as a version.
    m = re.search(r"Welcome to EQquasi\s+([0-9]+(?:\.[0-9]+)+)", out)
    if not m:
        raise RuntimeError(
            f"{exe} printed no recognisable version banner, so its build "
            f"cannot be checked against the source. Output was:\n{out[:2000]}")
    return m.group(1)


def check_binary_is_current(exe):
    """Refuse to run a benchmark against a binary older than the source.

    This lives here, in the path every benchmark and every reference
    generation goes through, rather than only in a test that reports the
    problem afterwards. A stale binary does not fail loudly: comparisons pass,
    because the oracle was frozen from the same stale build, and a reference
    blessed from it locks last week's physics in as ground truth. On
    2026-08-12 a full session of BP1002 step-over work ran on 1.7.0 while src
    was 1.10.0 -- across commits that changed the very multi-fault paths the
    case exercises -- and a half-generated reference had to be discarded.
    """
    declared = declared_version()
    if not declared:
        pytest.fail("src/globalvar.f90 declares no EQQUASI_VERSION, so the "
                    "binary cannot be checked against it.")
    try:
        built = binary_version(exe)
    except Exception as exc:
        pytest.fail(f"could not read {exe}'s version ({exc}). Refusing to run "
                    f"a benchmark against a binary of unknown provenance.")
    # The version string alone is not enough. A source edit that does not bump
    # the version -- which is most of them -- leaves src and bin agreeing on
    # 1.12.0 while the binary lacks the change: the check passes and the run
    # tests code that is no longer there. Compare build times too. This is the
    # same failure as the 1.7.0/1.10.0 skew, one step subtler.
    newest_src, newest_name = 0.0, ""
    src_dir = os.path.join(str(ROOT), "src")
    if os.path.isdir(src_dir):
        for f in os.listdir(src_dir):
            if f.endswith((".f90", ".F90")) or f == "makefile":
                m = os.path.getmtime(os.path.join(src_dir, f))
                if m > newest_src:
                    newest_src, newest_name = m, f
    if newest_src and os.path.getmtime(exe) < newest_src:
        import datetime as _dt
        fmt = lambda t: _dt.datetime.fromtimestamp(t).strftime("%Y-%m-%d %H:%M")
        pytest.fail(
            f"{exe} is older than src/{newest_name}: binary built "
            f"{fmt(os.path.getmtime(exe))}, source modified {fmt(newest_src)}. "
            f"Both report version {declared}, so the version check alone "
            f"cannot see this. Rebuild: bash install.eqquasi.sh -m <host>   "
            f"(utig, ubuntu, conda-linux, ...)")

    if declared != built:
        pytest.fail(
            f"bin/eqquasi is {built} but src/globalvar.f90 declares "
            f"{declared}. Running a benchmark against a stale binary tests "
            f"last week's code and would bless it as a reference. Rebuild: "
            f"bash install.eqquasi.sh -m <host>   (utig, ubuntu, conda-linux, ...)")


def run_case(compset, over, workdir):
    """Drive one case through the workflow and return its cycle0 directory.

    create.newcase -> case.setup -> run.sh, never the binary directly: run.sh
    owns the cycle loop, and bypassing it produces output shaped unlike any
    other run's (rule 19).
    """
    e = env()
    exe = installed_binary()
    if exe is None:
        # Not a skip. A skipped benchmark takes the entire e2e gate with it and
        # still reports green -- CI did exactly that: 2 passed, 18 skipped in
        # 8 seconds, having run no benchmark at all. Build first.
        pytest.fail(
            f"bin/eqquasi-{declared_version()} does not exist. Build before "
            f"running the e2e tier: bash install.eqquasi.sh -m <host>, which "
            f"installs the binary under its own version. There is no plain "
            f"bin/eqquasi by design.")

    check_binary_is_current(exe)

    def step(cmd, cwd):
        """Run a setup command, and on failure say what it printed.

        DEVNULL here made a CI failure undiagnosable: create.newcase exited 1
        and the only evidence was CalledProcessError. A test that hides the
        error it just caught is worse than no test.
        """
        r = subprocess.run(cmd, cwd=cwd, env=e, capture_output=True, text=True)
        if r.returncode != 0:
            pytest.fail(f"{' '.join(cmd)} failed (exit {r.returncode}) in {cwd}\n"
                        f"stdout:\n{r.stdout[-3000:]}\n"
                        f"stderr:\n{r.stderr[-3000:]}")

    # scratch/ is gitignored, so it does not exist in a fresh checkout and
    # create.newcase fails on the missing parent. Locally it is always there,
    # which is why this only ever failed in CI.
    os.makedirs(os.path.dirname(os.path.abspath(workdir)), exist_ok=True)
    shutil.rmtree(workdir, ignore_errors=True)
    # ignore_errors hides the one failure that matters. On NFS, a file still
    # open by a process -- ranks left behind when a previous run was killed --
    # is silly-renamed to .nfsXXXX rather than unlinked, so the directory
    # survives, create.newcase refuses it, and every category of every
    # remaining case fails on a comparison it never got to make. Twenty
    # failures, none of them saying "a stale process is holding the old run".
    if os.path.exists(workdir):
        leftovers = sorted(os.listdir(workdir))[:8]
        pytest.fail(
            f"{workdir} still exists after rmtree; create.newcase will refuse "
            f"it. Most likely a previous run's MPI ranks are still alive and "
            f"holding files open (NFS .nfs* entries). Check with "
            f"`ps -eo pid,cmd | grep eqquasi`, kill them, then remove the "
            f"directory. Leftovers: {leftovers}")
    # Absolute path, never PATH. The user's login environment appends an
    # OLDER EQquasi checkout's scripts/ to PATH (~/.bashrc EQQUASIROOT
    # pointing at 0.Dunyu/EQquasi), and a bare "create.newcase" resolved
    # there -- a copy that still expects case_input/ -- failing every e2e
    # case with FileNotFoundError while looking like a defect in THIS repo.
    # The gate must run this repo's tools, not whichever namesake PATH finds
    # first (rule 3: evaluated at the right configuration).
    step([str(ROOT / "scripts" / "create.newcase"), workdir, compset],
         str(ROOT))
    # Compsets that ship a geometry generator (the kink family) need it run
    # between create and setup: it writes input/bFault_Rough_Geometry.txt
    # from user_defined_params.py. Generic on the filename so the next
    # generator-bearing compset needs no harness change.
    import glob as _glob
    for gen in sorted(_glob.glob(os.path.join(workdir, "input",
                                              "generate*Geometry.py"))):
        step(["python3", gen, os.path.join(workdir, "input")], workdir)
    apply_overrides(os.path.join(workdir, "user_defined_params.py"), over)
    step(["./case.setup"], workdir)

    os.makedirs(os.path.join(workdir, "log"), exist_ok=True)
    with open(os.path.join(workdir, "log", "run.log"), "w") as log:
        proc = subprocess.Popen(["bash", "run.sh"], cwd=workdir, env=e,
                                stdout=log, stderr=subprocess.STDOUT,
                                start_new_session=True)
    t0 = time.time()
    while proc.poll() is None:
        if time.time() - t0 > RUN_TIMEOUT_S:
            proc.kill()
            pytest.fail(f"{compset} exceeded {RUN_TIMEOUT_S} s")
        time.sleep(10)
    if proc.returncode != 0:
        # Print the log, do not merely point at it. On a CI runner the file is
        # gone the moment the job ends, so "see workdir/run.log" is an
        # instruction nobody can follow: BP8 failed this way on GitHub for two
        # releases running and the only evidence available was the exit code.
        try:
            tail = open(os.path.join(workdir, "log", "run.log")).read()[-4000:]
        except OSError as exc:
            tail = f"(run.log unreadable: {exc})"
        pytest.fail(f"run.sh exited {proc.returncode} for {compset}.\n"
                    f"--- tail of {workdir}/log/run.log ---\n{tail}")

    # Cycles live under result/ since the case layout gained input/, result/,
    # scratch/, tool/ and log/. A run that failed leaves scratch/ and no
    # result/cycle0 at all, which is the point of the split.
    cyc = os.path.join(workdir, "result", "cycle0")
    assert os.path.isdir(cyc), (
        f"run.sh produced no result/cycle0/ for {compset}"
        + (f"; scratch/ is present with {len(os.listdir(os.path.join(workdir, 'scratch')))}"
           " files, so the cycle started and did not finish"
           if os.path.isdir(os.path.join(workdir, "scratch")) else ""))

    # A directory is not a result. When BP1002 hit STOP 508 in cycle 2, the
    # loop carried on and left two more cycleN/ directories behind, each with
    # one snapshot written at step 1 and no RUN SUMMARY -- output shaped like a
    # completed cycle, holding nothing. runInfo.json is written only after the
    # time loop finishes, so its presence, and a nonzero step count in it, is
    # what separates a run from a directory.
    info = os.path.join(cyc, "runInfo.json")
    if not os.path.exists(info):
        pytest.fail(
            f"{cyc} has no runInfo.json, so the time loop never completed for "
            f"{compset}. The run stopped early -- see {workdir}/run.log for "
            f"the reason (a STOP code, or the step cap).")
    try:
        steps = json.load(open(info)).get("steps_completed", 0)
    except (ValueError, OSError) as exc:
        pytest.fail(f"{info} is unreadable ({exc}); the run did not finish "
                    f"cleanly for {compset}.")
    if not steps:
        pytest.fail(f"{info} reports steps_completed = {steps} for {compset}: "
                    f"the run produced a cycle directory but took no steps.")
    return cyc


def read_any(path):
    """Numeric rows from a run (.txt/.dat, whitespace) or a reference (.csv)."""
    from seasio import read_array
    return read_array(path)


def station_files(d, prefix):
    return sorted(f for f in os.listdir(d)
                  if f.startswith(prefix) and f.endswith((".txt", ".dat", ".csv")))


def counterpart(ref_dir, run_dir, name):
    """The run file matching a reference file, whatever its extension."""
    stem = os.path.splitext(name)[0]
    # .nc must be here: a fault.*.nc reference otherwise has no counterpart and
    # every snapshot reports "not produced by the run".
    for ext in (".nc", ".dat", ".txt", ".csv"):
        p = os.path.join(run_dir, stem + ext)
        if os.path.exists(p):
            return p
    return None


def load_summary(ref_dir):
    p = os.path.join(ref_dir, "summary.json")
    return json.load(open(p)) if os.path.exists(p) else None


# --------------------------------------------------------------------------
# Categories: what kinds of file a reference can hold, how each is compared,
# and how a result is reported. One method per category, chosen by what the
# file *is*, not by which benchmark produced it.

CATEGORIES = {
    "snapshot": (("fault.0*.nc",),                   "netcdf"),
    "restart":  (("fault.r.nc",),                    "netcdf"),
    # The last-step whole-model displacement. Stored as float64 (the solver's
    # dp) since 1.24.8; compared, like fault.*.nc, at compare_netcdf's 1e-10.
    "disp":     (("disp.0*.nc",),                    "netcdf"),
    "onfault":  (("fltst_strk*",),                   "series"),
    "offfault": (("srfst_strk*",),                   "series"),
    "global":   (("global.dat", "global.csv"),       "series"),
    # Row 6 (2026-09-25): time + one peak-V column per fault. Additive to
    # global.dat (whose column 2 stays the all-fault max); only the
    # references re-blessed additively since then hold it.
    "perfault": (("peak_sliprate_per_fault.dat",),    "series"),
    "profile":  (("*_strike.dat", "*_strike.csv",
                  "*_depth.dat",  "*_depth.csv"),    "series"),
    "cplot":    (("cplot_EQquasi.*", "cplot_ft2_EQquasi.*",
                  "cplot_ruptarea_trac_slip.*",
                  "cplot_ft2_ruptarea_trac_slip.*"), "series"),
    "scalars":  (("summary.json",),                  "summary"),
}


def category_files(ref_dir, category):
    """Sorted reference files in one category, [] when the case has none."""
    import glob as _g
    pats, _ = CATEGORIES[category]
    out = []
    for pat in pats:
        out += _g.glob(os.path.join(ref_dir, pat))
    return sorted(set(out))


def manifest(ref_dir):
    """{category: [files]} for everything this reference actually holds.

    A case is defined by its files, so the manifest is the case: BP8 has
    profiles and no off-fault stations, the step-over has two fault slabs, BP5
    has a cplot with real rupture times. Nothing here names a benchmark.
    """
    return {c: f for c in CATEGORIES
            for f in (category_files(ref_dir, c),) if f}


def _scaled_tol(ref, rtol, floor_rtol, scale):
    """Elementwise tolerance: the entry's own size, plus a floor from `scale`.

        rtol * |ref|          -- the entry itself
        floor_rtol * scale    -- the size of the QUANTITY it belongs to

    `scale` must come from a sibling of the same physical kind -- the other
    components of the same vector -- and never from the file as a whole. An
    earlier version used the file's maximum, which in a station file is time
    in seconds (2.6e6); every slip-rate entry then sat below the resulting
    threshold and was dropped from comparison entirely, so halving BP8's slip
    rate at all 5302 rows passed with "worst diff = 0.0e+00". Nothing is
    exempted now; a component that vanishes by symmetry is merely judged
    against the component that does not.

    ATOL is the machine-noise floor for the case the sibling rule cannot
    reach: a quantity group that is zero to machine precision in ref AND run.
    First hit by stepover's station at the locked VW centre -- after 101
    steps its slips and rates are +-1e-27 MPI-reduction noise on both sides,
    the group scale is itself noise, and rtol on noise fails a third of the
    file. 1e-18 sits ~9 orders above that noise and ~9 below the smallest
    physical signal any comparison carries (creep, 1e-9 m/s). It cannot mask
    the BP8 failure mode above: a halved 1e-9-scale signal differs by ~5e-10,
    twenty-eight orders past this floor.
    """
    import numpy as _np
    return rtol * _np.abs(ref) + floor_rtol * scale + ATOL_MACHINE_NOISE


# Module constant, not a per-call default: every caller should mean the same
# thing by "zero to machine precision".
ATOL_MACHINE_NOISE = 1e-18


def compare_arrays(ref, run, rtol=1e-9, floor_rtol=1e-9, scale=None):
    """(ok, n_bad, worst_index, worst_scaled) for two same-shaped arrays.

    One rule for tables and netCDF alike: they used to disagree, tables
    comparing against a fixed 1e-12 floor and netCDF demanding bit equality,
    so the same floating-point noise was a pass in one and a failure in the
    other.
    """
    import numpy as _np
    ref = _np.asarray(ref, float)
    run = _np.asarray(run, float)
    # A NaN compares False against every tolerance, so `diff > tol` silently
    # accepts it. Catch non-finite values explicitly and on both sides: a run
    # that returns NaN is the single most important thing this gate can say.
    bad_finite = _np.isfinite(ref) != _np.isfinite(run)
    bad_finite |= ~_np.isfinite(run) & ~_np.isfinite(ref) & (
        _np.nan_to_num(run, nan=0.0, posinf=1.0, neginf=-1.0)
        != _np.nan_to_num(ref, nan=0.0, posinf=1.0, neginf=-1.0))
    if scale is None:
        scale = float(_np.max(_np.abs(ref[_np.isfinite(ref)]))) if ref.size else 0.0
    diff = _np.abs(run - ref)
    tol = _scaled_tol(ref, rtol, floor_rtol, scale)
    bad = (_np.isfinite(diff) & (diff > tol)) | bad_finite
    ratio = _np.where(_np.isfinite(diff), diff / _np.maximum(tol, 1e-300), _np.inf)
    k = int(_np.argmax(_np.where(_np.isfinite(ratio), ratio, _np.finfo(float).max)))
    idx = _np.unravel_index(k, ref.shape)
    denom = float(_np.max(_np.abs(ref[_np.isfinite(ref)]))) if ref.size else 0.0
    worst = float(_np.max(diff[_np.isfinite(diff)]) / (denom or 1.0)) if diff.size else 0.0
    return (not bad.any()), int(bad.sum()), idx, worst


def column_group_scales(ref, header):
    """Per-column scale: the largest magnitude among that column's siblings.

    SEAS station columns name a vector's components with a trailing _2 / _3
    (slip_2/slip_3, slip_rate_2/slip_rate_3, shear_stress_2/shear_stress_3).
    Stripping that suffix groups them, so BP8's dip components -- which are
    zero by symmetry and carry only rounding -- are measured against their
    strike siblings instead of against themselves. Without a header every
    column stands alone, which is right for cplot and global.dat, whose
    near-zero entries sit inside columns that do carry signal.
    """
    import numpy as _np
    col_max = _np.max(_np.abs(_np.where(_np.isfinite(ref), ref, 0.0)), axis=0)
    if not header or len(header) != ref.shape[1]:
        return col_max
    groups = {}
    for j, name in enumerate(header):
        g = name[:-2] if name.endswith(("_2", "_3")) else name
        groups.setdefault(g, []).append(j)
    out = col_max.copy()
    for cols in groups.values():
        out[cols] = col_max[cols].max()
    return out


# Columns the SEAS output format stores as log10, by header name. They are
# compared in linear space: BP8's slip_rate_3 is the dip slip rate on a fault
# with no dip motion, so it wanders between 1e-30 and 1e-20 m/s -- all of it
# zero next to the 1e-12 m/s strike rate -- but in log10 that reads as a
# difference of 4, and no tolerance on the logarithm can tell that from a real
# change. Taking 10**x first puts the comparison back where zero means zero.
LOG10_COLUMNS = ("slip_rate_2", "slip_rate_3")


def _header_of(*paths, ncols=None):
    """Column names from whichever path is a .csv with a usable header."""
    for p in paths:
        if not isinstance(p, str) or not p.endswith(".csv"):
            continue
        try:
            with open(p) as fh:
                names = fh.readline().strip().split(",")
        except OSError:
            continue
        if ncols is not None and len(names) != ncols:
            continue
        return names
    return None


def sibling_component_scale(path, ncols):
    """Largest magnitude across the sibling files of a vector component.

    BP8's section 4.3 profiles put each component in its own file --
    shear_stress_2_depth.csv holds strike, shear_stress_3_depth.csv holds dip.
    A file that is entirely one component has no in-file sibling to be scaled
    against, and dip shear on this fault is 1e-16 Pa of rounding, so judged
    against itself it fails every time while carrying no information. Look up
    the other component's file and use its magnitude instead.

    Returns None when the file is not a component file or no sibling exists,
    in which case the caller falls back to per-column scaling.
    """
    import glob as _g
    import numpy as _np
    import os as _os
    import re as _re
    base = _os.path.basename(path)
    m = _re.match(r"^(.*?)_([23])_(.*)$", base)
    if not m:
        return None
    stem, _comp, tail = m.groups()
    best = 0.0
    for sib in _g.glob(_os.path.join(_os.path.dirname(path) or ".",
                                     f"{stem}_[23]_{tail}")):
        try:
            a = _np.atleast_2d(_np.asarray(read_any(sib), float))
        except Exception:
            continue
        if a.shape[1] != ncols:
            continue
        finite = a[_np.isfinite(a)]
        if finite.size:
            best = max(best, float(_np.max(_np.abs(finite))))
    return best or None


def column_printed_rtol(path, ncols, rtol):
    """Per-column relative tolerance, floored at the file's own printed precision.

    A text table can only resolve differences it has the digits to express.
    The station files are written `E21.13,6E15.7`: seven significant digits in
    every column but the first. Two runs that differ by MPI reduction ordering
    (~2e-14 relative) print identical text almost everywhere, but a value
    sitting on a rounding boundary flips its last digit -- 1.383710E-05 against
    1.383711E-05, one entry in 31381. Demanding rtol=1e-9 of a column carrying
    1e-7 of precision asks the file for digits it never wrote, so the gate
    fails at random on reruns rather than on regressions.

    So the floor is read off the file: `d` significant digits give one unit in
    the last place of 10**-(d-1) relative, doubled for the case where a value
    rounds away from its neighbour rather than toward it. Nothing is loosened
    beyond what the text cannot say -- the netCDF files carry full double
    precision and are still compared at 1e-9, and a regression too small to
    change the seventh digit of a station file is not one this tier could ever
    have detected.

    Only exponential-notation tokens are measured. A plain "1.0" would read as
    two significant digits and floor the column at 1e-1, which would be a real
    loosening; columns without exponential tokens keep the caller's rtol.
    """
    import numpy as _np
    import re as _re
    tok = _re.compile(r"^[+-]?(\d*)\.?(\d*)[EeDd][+-]?\d+$")
    digits = [None] * ncols
    try:
        with open(path) as fh:
            for line in fh:
                parts = line.split()
                if len(parts) != ncols:
                    continue
                for j, p in enumerate(parts):
                    m = tok.match(p)
                    if not m:
                        continue
                    sig = (m.group(1) + m.group(2)).lstrip("0")
                    n = len(sig)
                    if n and (digits[j] is None or n < digits[j]):
                        digits[j] = n          # least precise token wins
    except OSError:
        return _np.full(ncols, float(rtol))

    out = _np.full(ncols, float(rtol))
    for j, d in enumerate(digits):
        if d:
            out[j] = max(out[j], min(2.0 * 10.0 ** -(d - 1), 1e-5))
    return out


# Round-off floors of the fltst_ out-of-plane columns (owner decision, board
# row 1, 2026-11-18/19), in the un-logged units of each column: slip_3 in m,
# slip_rate_3 in m/s (the 11-column file holds its log10), shear_stress_3 in
# MPa. An entry whose un-logged values differ by no more than the floor is
# skipped; so is a log10 entry whose un-logged values both sit under it.
STATION_FLOORS = {"slip_3": 1e-15, "slip_rate_3": 1e-12, "shear_stress_3": 1e-9}
# (0-based column, floor name, holds log10) by fltst_ column count
# (src/library_output.f90): 9 for BP5/BP7, 11 for BP8.
STATION_FLOOR_COLUMNS = {
    11: ((2, "slip_3", False), (4, "slip_rate_3", True), (6, "shear_stress_3", False)),
    9: ((2, "slip_3", False), (4, "slip_rate_3", False), (6, "shear_stress_3", False)),
}


def station_floor_exempt(gpath, g, r):
    """Boolean mask of the entries of two fltst_ tables that sit inside a
    round-off floor. Only a file named fltst_* in a 9- or 11-column layout is
    touched (no srfst_ floor); every other file gets an all-False mask. An
    entry is exempt only when both its reference and run raw values (the
    log10 values, un-exponentiated) are finite; g and r must be raw."""
    import numpy as _np
    exempt = _np.zeros(g.shape, bool)
    cols = STATION_FLOOR_COLUMNS.get(g.shape[1])
    if cols is None or not os.path.basename(gpath).startswith("fltst_"):
        return exempt
    for j, name, log10 in cols:
        lin = (lambda a: 10.0 ** a) if log10 else (lambda a: a)
        with _np.errstate(invalid="ignore", over="ignore"):
            exempt[:, j] = (_np.isfinite(g[:, j]) & _np.isfinite(r[:, j])
                            & (_np.abs(lin(r[:, j]) - lin(g[:, j])) <= STATION_FLOORS[name]))
    return exempt


def compare_series(gpath, rpath, rtol=1e-9, floor_rtol=1e-9):
    """(ok, message) for two numeric tables. Full array, not sampled.

    Relative, with a per-quantity floor. A fixed absolute tolerance cannot
    work here: one file holds time in seconds (1e6), moment rate (1e11), slip
    in metres (1e-2) and log10 slip rates near -30.
    """
    import numpy as _np
    g, r = read_any(gpath), read_any(rpath)
    if g.shape != r.shape:
        return False, f"shape {r.shape} vs reference {g.shape}"
    g, r = _np.atleast_2d(_np.array(g, float)), _np.atleast_2d(_np.array(r, float))

    g_raw, r_raw = g.copy(), r.copy()
    header = _header_of(gpath, rpath, ncols=g.shape[1])
    if header:
        for j, name in enumerate(header):
            if name in LOG10_COLUMNS:
                g[:, j] = 10.0 ** g[:, j]
                r[:, j] = 10.0 ** r[:, j]

    scale = column_group_scales(g, header)
    sib = sibling_component_scale(gpath, g.shape[1])
    if sib is not None:
        scale = _np.maximum(scale, sib)
    # A log10 column was exponentiated above, so its printed precision no
    # longer describes the numbers being compared; measure before that and it
    # would be wrong, so those columns simply keep the caller's rtol.
    col_rtol = column_printed_rtol(gpath, g.shape[1], rtol)
    if header:
        for j, name in enumerate(header):
            if name in LOG10_COLUMNS:
                col_rtol[j] = rtol

    # An entry inside a round-off floor is made equal to its reference, so
    # compare_arrays skips it; the scale was measured before this.
    exempt = station_floor_exempt(gpath, g_raw, r_raw)
    r = _np.where(exempt, g, r)

    ok, nbad, idx, worst = compare_arrays(g, r, col_rtol, floor_rtol, scale=scale)
    if not ok:
        i, j = idx
        return False, (f"{nbad} of {g.size} entries outside "
                       f"rtol={col_rtol[j]:g}; "
                       f"worst at row {i}, column {j+1}: "
                       f"{r[i, j]:.8g} vs {g[i, j]:.8g}")
    return True, (f"{g.shape[0]} rows x {g.shape[1]} cols, "
                  f"worst diff = {worst:.1e} of the file's scale")


def variable_group(name):
    """The family a fault variable belongs to: its vector, not itself.

    A component that vanishes by symmetry has to be judged against the vector
    it is a component of. `slipn` is fault opening, forbidden by the contact
    condition, so it sits at 4.6e-26 m while `slips` carries 3.7e-2 m of real
    slip; `vzm` is 3.5e-35 m/s beside `vxm`'s 5e-10.
    """
    for suffix in ("_strike", "_dip", "_normal"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    if len(name) > 1 and name[:-1] in ("slip",) and name[-1] in "sdn":
        return name[:-1]
    if len(name) == 3 and name[0] == "v" and name[1] in "xyz" and name[2] in "ms":
        return f"v_{name[2]}"
    return name


def compare_netcdf(gpath, rpath, rtol=1e-10, floor_rtol=1e-10):
    """(ok, message) per fault, so a missing or aliased fault is named.

    Same tolerance rule as compare_series, with each variable's vector as the
    scale. Bit equality is not available: the same case, same binary, same
    host, run twice on 3 MPI ranks differs by ~2e-14 from reduction ordering
    alone.
    """
    import numpy as _np
    import netCDF4 as _nc
    g, r = _nc.Dataset(gpath), _nc.Dataset(rpath)
    gd = {k: len(v) for k, v in g.dimensions.items()}
    rd = {k: len(v) for k, v in r.dimensions.items()}
    if gd != rd:
        return False, f"dimensions {rd} vs reference {gd}"

    # nid_* are the per-fault index maps -- the very thing bp1002 exists to
    # check -- so they are compared exactly rather than skipped.
    names = [v for v in g.variables]
    for v in names:
        if v.startswith("nid_"):
            a = _np.asarray(g[v][:]); b = _np.asarray(r[v][:])
            if a.shape != b.shape or not _np.array_equal(a, b):
                return False, f"index map {v} differs from the reference"

    names = [v for v in names if not v.startswith("nid_")]
    group_scale = {}
    for v in names:
        grp = variable_group(v)
        a = _np.asarray(g[v][:], float)
        m = float(_np.max(_np.abs(a[_np.isfinite(a)]))) if a.size else 0.0
        group_scale[grp] = max(group_scale.get(grp, 0.0), m)

    worst, where, nbad = 0.0, "", 0
    for v in names:
        a, b = _np.asarray(g[v][:], float), _np.asarray(r[v][:], float)
        nf = a.shape[0] if a.ndim == 3 else 1
        for ift in range(nf):
            x = a[ift] if a.ndim == 3 else a
            y = b[ift] if b.ndim == 3 else b
            ok, n, _idx, rel = compare_arrays(
                x, y, rtol, floor_rtol, scale=group_scale[variable_group(v)])
            nbad += n
            if rel > worst:
                worst = rel
                where = f"{v}{f'[fault {ift}]' if nf > 1 else ''}"
    if nbad:
        return False, (f"{nbad} entries outside rtol={rtol:g}; worst is "
                       f"{worst:.3e} of the field's scale in {where}")
    return True, (f"{len(gd)} dims, {len(g.variables)} variables, "
                  f"worst diff = {worst:.1e} of the field's scale")


COMPARATORS = {"series": compare_series, "netcdf": compare_netcdf}


# ---------------------------------------------------------------------------
# Full-cycle physical metrics (board PATHWAY_FORWARD.md, D6, owner option A,
# 2026-11-24). Definitions and tolerances are the board's; do not loosen.
#
# Event: contiguous samples of global.dat column 2 (peak slip rate) at or above
# EVENT_THRESHOLD. Onset and end: log-linear threshold crossing. Moment:
# trapezoid integral of column 3 over the event window. Coseismic slip: the
# fltst station with the largest slip, hypot(columns 2, 3) at end minus onset.
# A window touching the first or last sample is truncated: its onset or end is
# NaN and is never compared (D6 windows).
EVENT_THRESHOLD = 1e-3  # m/s
D6_TOLERANCES = {   # bp5 / bp7, solver 4 against solver 3 (round-off class)
    "onset_s": 4e-15,   # absolute, seconds
    "end_s": 1e-10,     # absolute, seconds
    "peak": 6e-14,      # relative
    "moment": 2e-14,    # relative
    "slip": 2e-14,      # relative
}


def load_columns(path):
    """Whitespace-separated float table; empty or ragged files raise."""
    rows = [[float(x) for x in ln.split()] for ln in open(path) if ln.split()]
    if not rows:
        raise ValueError(f"empty table {path}")
    if len({len(r) for r in rows}) != 1:
        raise ValueError(f"ragged table {path}")
    return rows


def log_crossing(t0, v0, t1, v1, thr):
    """Time at which log(v) (linear in t between the samples) equals log(thr)."""
    if v0 <= 0 or v1 <= 0:
        raise ValueError(f"non-positive slip rate at a crossing: {v0}, {v1}")
    return t0 + (math.log(thr) - math.log(v0)) / (math.log(v1) - math.log(v0)) * (t1 - t0)


def _lin(t, ts, ys, i):
    return ys[i] + (t - ts[i]) / (ts[i + 1] - ts[i]) * (ys[i + 1] - ys[i])


def find_events(g, thr=EVENT_THRESHOLD):
    """Events of a global.dat table (rows: time, peak slip rate, moment, ...).

    Returns dicts with onset, end (NaN when truncated), t0, t1 (the window used),
    peak and moment.
    """
    t = [r[0] for r in g]
    v = [r[1] for r in g]
    m = [r[2] for r in g]
    n = len(t)
    out = []
    i = 0
    while i < n:
        if v[i] >= thr:
            a = i
            while i + 1 < n and v[i + 1] >= thr:
                i += 1
            b = i
            ts, vs = a == 0, b == n - 1
            on = math.nan if ts else log_crossing(t[a - 1], v[a - 1], t[a], v[a], thr)
            en = math.nan if vs else log_crossing(t[b], v[b], t[b + 1], v[b + 1], thr)
            tw = ([] if ts else [on]) + t[a:b + 1] + ([] if vs else [en])
            mw = (([] if ts else [_lin(on, t, m, a - 1)]) + m[a:b + 1]
                  + ([] if vs else [_lin(en, t, m, b)]))
            mom = sum(0.5 * (mw[k] + mw[k + 1]) * (tw[k + 1] - tw[k]) for k in range(len(tw) - 1))
            out.append(dict(onset=on, end=en, t0=tw[0], t1=tw[-1],
                            peak=max(v[a:b + 1]), moment=mom))
        i += 1
    return out


def station_slip_over(fltst_path, t0, t1):
    """hypot(fltst cols 2, 3) at t1 minus t0, interpolated; NaN if outside the series."""
    f = load_columns(fltst_path)
    ts = [r[0] for r in f]
    s = [math.hypot(r[1], r[2]) for r in f]
    if t0 < ts[0] or t1 > ts[-1]:
        return math.nan

    def at(t):
        return _lin(t, ts, s, min(max(bisect.bisect_right(ts, t) - 1, 0), len(ts) - 2))
    return at(t1) - at(t0)


def full_cycle_metrics(result_dir, thr=EVENT_THRESHOLD):
    """Per-event metrics of a result/cycle0 directory (global.dat + fltst_*.txt).

    Each event also carries `station` (largest coseismic slip in THIS run) and
    `slip`. NaN-slip stations (truncated windows) never win.
    """
    result_dir = str(result_dir)
    events = find_events(load_columns(os.path.join(result_dir, "global.dat")), thr)
    stations = sorted(glob.glob(os.path.join(result_dir, "fltst_*.txt")))
    if not stations:
        raise FileNotFoundError(f"no fltst_*.txt in {result_dir}")
    for ev in events:
        slips = {p: station_slip_over(p, ev["t0"], ev["t1"]) for p in stations}
        ok = {p: x for p, x in slips.items() if not math.isnan(x)}
        best = max(ok, key=ok.get) if ok else None
        ev["station"] = os.path.basename(best) if best else None
        ev["slip"] = ok[best] if best else math.nan
    return events


def _rel(a, b):
    return 0.0 if a == b else abs(b - a) / max(abs(a), abs(b))


def compare_full_cycle_metrics(ref, run, tol=D6_TOLERANCES):
    """(ok, message) comparing two `full_cycle_metrics` lists at the D6 tolerances.

    NaN-safe: a failure is `not (err <= bound)`. Both-NaN (truncated) onset, end
    or slip is skipped; NaN on one side only fails.
    """
    if len(ref) != len(run):
        return False, f"event count {len(run)} != reference {len(ref)}"
    bad = []
    for k, (a, b) in enumerate(zip(ref, run)):
        pairs = (("onset_s", abs(b["onset"] - a["onset"]), a["onset"], b["onset"]),
                 ("end_s", abs(b["end"] - a["end"]), a["end"], b["end"]),
                 ("peak", _rel(a["peak"], b["peak"]), a["peak"], b["peak"]),
                 ("moment", _rel(a["moment"], b["moment"]), a["moment"], b["moment"]),
                 ("slip", _rel(a["slip"], b["slip"]), a["slip"], b["slip"]))
        for key, err, x, y in pairs:
            if math.isnan(x) and math.isnan(y):
                continue
            if not (err <= tol[key]):
                bad.append(f"event {k} {key}: {err:.3e} > {tol[key]:g}")
        if a["station"] != b["station"]:
            bad.append(f"event {k} station {b['station']} != reference {a['station']}")
    if bad:
        return False, "; ".join(bad)
    return True, f"{len(ref)} event(s) within D6 tolerances"
