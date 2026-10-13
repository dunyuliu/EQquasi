"""End-to-end: run each case, compare every reference file, SUCCESS or FAIL.

Replaces five per-benchmark files (~950 lines) that differed only in a compset
name, a few parameter overrides and which reference to diff against -- and had
each grown its own comparison logic, so BP8 got full-array diffs and BP5 got
none, for no reason beyond which file was edited last.

A case is defined by the files its reference holds. `cases.manifest()` sorts
them into categories, and each category has one comparison method:

    snapshot  fault.0*.nc        netcdf, per fault, relative to each variable
    restart   fault.r.nc         netcdf
    onfault   fltst_strk*        full series diff
    offfault  srfst_strk*        full series diff
    global    global.dat/.csv    full series diff
    profile   *_strike/_depth    full series diff  (BP8 only writes these)
    cplot     cplot_EQquasi.*    full series diff
    scalars   summary.json       the named values that reference calls important

Nothing here names a benchmark to decide what to check, so a new benchmark is a
row in `cases.CASES` plus a reference directory.

    pytest -m e2e_ci tests/      what CI runs on every push: two 101-step
                                 smokes (test.bp5.qdc.2000,
                                 test.stepover.qdc.1000 for rule 12)
    pytest -m e2e_fast tests/    the periodic/local sweep: every fast case,
                                 101-step and BP8's 8000-step alike, plus the
                                 clean build
    pytest -m e2e tests/         adds the full tiers (a full BP5 cycle, etc.)
"""

import json
import os
import shutil
import subprocess

import pytest

from conftest import ROOT
import cases as C

pytestmark = pytest.mark.e2e

WORK_ROOT = os.path.join(str(ROOT), "scratch")

# The two fast cases CI itself runs on every push (owner decision 2026-09-23,
# PATHWAY_FORWARD.md row 9): a plain BP5 smoke, and the project's only
# ntotft > 1 case in the gate (rule 12). Everything else under e2e_fast --
# including BP8's ~8000-step case, ~13 min alone -- moves to the
# periodic/local sweep and is never selected by CI.
CI_FAST_NAMES = {"test.bp5.qdc.2000", "test.stepover.qdc.1000"}


def _marks(c):
    if c[4] != "fast":
        return []
    marks = [pytest.mark.e2e_fast]
    if c[0] in CI_FAST_NAMES:
        marks.append(pytest.mark.e2e_ci)
    return marks


ALL = [pytest.param(c, marks=_marks(c)) for c in C.CASES]
IDS = [C.case_id(c) for c in C.CASES]


@pytest.fixture(scope="session")
def runs():
    """One run per case, shared across the checks that read it."""
    return {}


def _run(runs, case):
    cid = C.case_id(case)
    if cid not in runs:
        name, compset, over, ref, _ = case
        ref_dir = C.reference_dir(name, ref)
        if not os.path.isdir(ref_dir):
            pytest.skip(f"no reference at {ref_dir}")
        runs[cid] = (C.run_case(compset, over,
                                os.path.join(WORK_ROOT, f"e2e.{cid}")), ref_dir)
    return runs[cid]


@pytest.mark.parametrize("case", ALL, ids=IDS)
@pytest.mark.parametrize("category", [c for c in C.CATEGORIES if c != "scalars"])
def test_category_matches_reference(runs, case, category):
    """Compare every reference file in one category, and report each by name.

    Reported per file rather than as a single verdict: "bp5 failed" sends you
    looking, "fltst_strk000dp010.txt: max|diff| = 3.1e-04 at column 6" does not.
    """
    run_dir, ref_dir = _run(runs, case)
    files = C.category_files(ref_dir, category)
    if not files:
        pytest.skip(f"this reference holds no {category} files")

    compare = C.COMPARATORS[C.CATEGORIES[category][1]]
    lines, failed = [], False
    for g in files:
        name = os.path.basename(g)
        r = C.counterpart(ref_dir, run_dir, name)
        if r is None:
            lines.append(f"FAIL    {name}: not produced by the run")
            failed = True
            continue
        try:
            ok, msg = compare(g, r)
        except Exception as exc:                      # a comparator that dies
            ok, msg = False, f"{type(exc).__name__}: {exc}"
        lines.append(f"{'SUCCESS' if ok else 'FAIL   '} {name}: {msg}")
        failed |= not ok

    report = "\n".join(lines)
    print(f"\n{C.case_id(case)} / {category}\n{report}")
    assert not failed, f"\n{C.case_id(case)} / {category}\n{report}"


@pytest.mark.parametrize("case", ALL, ids=IDS)
def test_summary_scalars_match(runs, case):
    """summary.json is where a reference states what it considers important."""
    run_dir, ref_dir = _run(runs, case)
    s = C.load_summary(ref_dir)
    if not s:
        pytest.skip("this reference has no summary.json")

    lines, failed = [], False

    def check(label, got, want, rel=1e-4):
        nonlocal failed
        # Absolute floor: peak slip rate is legitimately ~0 before nucleation,
        # and a pure relative test on zero can never pass.
        ok = abs(got - want) <= abs(want) * rel + 1e-12
        lines.append(f"{'SUCCESS' if ok else 'FAIL   '} {label}: "
                     f"{got:.6g} vs {want:.6g}")
        failed |= not ok

    if "run" in s and "steps" in s["run"]:
        info = json.load(open(os.path.join(run_dir, "runInfo.json")))
        check("steps", info["steps_completed"], s["run"]["steps"], rel=0.02)
    if "event" in s and "peak_Vmax_m_s" in s["event"]:
        # global.dat column 2 is LINEAR m/s for bp != 8; only the BP8 writer
        # applies dlog10. Do not exponentiate it.
        g = C.read_any(os.path.join(run_dir, "global.dat"))
        check("peak slip rate", float(g[:, 1].max()),
              s["event"]["peak_Vmax_m_s"], rel=1e-3)

    for st, vals in s.items():
        if not st.startswith(("+", "-")) or not isinstance(vals, dict):
            continue
        r = C.counterpart(ref_dir, run_dir, f"fltst_strk{st}.dat")
        if not r:
            continue
        a = C.read_any(r)
        if "slip_mm" in vals:
            check(f"{st} slip (mm)", a[-1, 1] * 1000, vals["slip_mm"])
        if "p_late_MPa" in vals and a.shape[1] > 7:
            check(f"{st} late p (MPa)", a[-1, 7], vals["p_late_MPa"])

    if not lines:
        pytest.skip("summary.json holds nothing this checks")
    report = "\n".join(lines)
    print(f"\n{C.case_id(case)} / scalars\n{report}")
    assert not failed, f"\n{C.case_id(case)} / scalars\n{report}"


@pytest.mark.parametrize("case", ALL, ids=IDS)
def test_output_is_uploadable(runs, case):
    """A result must also pass the submission validator, not merely match.

    Only meaningful where the section 4.3 profiles exist, which is BP8.
    """
    run_dir, ref_dir = _run(runs, case)
    if not C.category_files(ref_dir, "profile"):
        pytest.skip("no section 4.3 profiles in this reference")

    out = os.path.join(WORK_ROOT, f"e2e.sub.{C.case_id(case)}")
    subprocess.run(["python3", os.path.join(str(ROOT), "scripts",
                                            "resampleBP8Profiles.py"),
                    run_dir, out], check=True, stdout=subprocess.DEVNULL)
    v = subprocess.run(["python3", os.path.join(str(ROOT), "scripts",
                                                "checkBP8Submission"), out],
                       capture_output=True, text=True)
    print("\n" + v.stdout[-1500:])
    assert v.returncode == 0, "checkBP8Submission rejected the run"


@pytest.mark.e2e_fast
def test_restart_rejects_mismatched_dip_spacing():
    """PATHWAY_FORWARD.md row 5: a coupled restart's dip spacing must equal dx.

    EQdyna writes fault.dyna.r.nc (copied to fault.r.nc by run.sh) in
    EQquasi's multi-fault restart layout, indexing the dip axis with its own
    dz; EQquasi's own mesh never has an independent dz (src/mesh4num.f90,
    src/meshgen.f90, src/func_lib.f90 all set dz = dx by construction), so
    netcdf_read_on_fault_restart (src/netcdf_io.f90) indexes the same axis
    with dx. A dz != dx coupled restart would silently misalign the fault
    grid along dip: same (nfxMax, nfzMax, ntotft) shape either way, so no
    existing shape check catches it. netcdf_write_on_fault now stamps a
    "dip_node_spacing" global attribute and netcdf_read_on_fault_restart
    checks it before using infile2 (PROJECT_RULES.md rule 2).

    This exercises the real Fortran code path (rule 3): a fresh istart=2
    case, with result/cycle0 pre-seeded from a genuine cycle-0 run's
    fault.r.nc -- untouched for the pass case, with the attribute edited for
    the fail case -- then `bash run.sh` runs cycle 2 for real, which reads
    scratch/fault.r.nc via netcdf_read_on_fault_restart exactly as a true
    second cycle would (src/eqquasi.f90's checkAndReport, icstart > 1).
    Pre-seeding result/cycle0 before the case's first run.sh call is not the
    "inputs frozen after a cycle exists" door in scripts/case.setup: that
    guards against re-running case.setup, not against what run.sh's own
    cycle loop finds on disk.
    """
    netCDF4 = pytest.importorskip("netCDF4")
    compset = "test.stepover.qdc.1000"

    seed_dir = C.run_case(compset, {"istart": 1, "iend": 1},
                          os.path.join(WORK_ROOT, "e2e.restart_guard_seed"))
    seed_fault_r = os.path.join(seed_dir, "fault.r.nc")
    assert os.path.isfile(seed_fault_r), f"no fault.r.nc produced at {seed_dir}"

    with netCDF4.Dataset(seed_fault_r) as d:
        assert "dip_node_spacing" in d.ncattrs(), (
            f"{seed_fault_r} carries no dip_node_spacing attribute; "
            "netcdf_write_on_fault should stamp it (src/netcdf_io.f90)")
        seeded_spacing = float(d.getncattr("dip_node_spacing"))

    def seed_and_run(name, mutate):
        workdir = os.path.join(WORK_ROOT, f"e2e.{name}")
        e = C.env()
        shutil.rmtree(workdir, ignore_errors=True)
        subprocess.run([str(ROOT / "scripts" / "create.newcase"), workdir, compset],
                       cwd=str(ROOT), env=e, check=True, capture_output=True)
        C.apply_overrides(os.path.join(workdir, "user_defined_params.py"),
                          {"istart": 2, "iend": 2})
        subprocess.run(["./case.setup"], cwd=workdir, env=e, check=True,
                       capture_output=True)
        seeded_cycle0 = os.path.join(workdir, "result", "cycle0")
        os.makedirs(seeded_cycle0, exist_ok=True)
        seeded_fault_r = os.path.join(seeded_cycle0, "fault.r.nc")
        shutil.copy(seed_fault_r, seeded_fault_r)
        mutate(seeded_fault_r)
        os.makedirs(os.path.join(workdir, "log"), exist_ok=True)
        log_path = os.path.join(workdir, "log", "run.log")
        with open(log_path, "w") as log:
            r = subprocess.run(["bash", "run.sh"], cwd=workdir, env=e,
                               stdout=log, stderr=subprocess.STDOUT)
        return r.returncode, open(log_path).read()

    # (a) the genuine, unmutated restart file: dz == dx by construction, so
    # this must still succeed exactly as it did before the guard existed.
    rc_ok, log_ok = seed_and_run("restart_guard_ok", lambda p: None)
    assert rc_ok == 0, (
        f"a genuine, unmutated restart file was rejected: exit {rc_ok}\n"
        f"{log_ok[-3000:]}")

    # (b) a mismatched dip_node_spacing: what a coupled EQdyna restart with
    # dz != dx looks like once EQdyna stamps the same attribute.
    def corrupt(path):
        with netCDF4.Dataset(path, "r+") as d:
            d.setncattr("dip_node_spacing", seeded_spacing * 1.5)

    rc_bad, log_bad = seed_and_run("restart_guard_mismatch", corrupt)
    assert rc_bad != 0, "a mismatched dip-axis spacing did not fail loudly"
    assert "dip-axis node spacing" in log_bad and "does not match" in log_bad, (
        f"failed, but not with the expected message:\n{log_bad[-3000:]}")
