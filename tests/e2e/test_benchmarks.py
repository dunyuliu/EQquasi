"""End-to-end: run each case, compare every reference file, SUCCESS or FAIL.

Replaces five per-benchmark files (~950 lines) that differed only in a compset
name, a few parameter overrides and which reference to diff against -- and had
each grown its own comparison logic, so BP8 got full-array diffs and BP5 got
none, for no reason beyond which file was edited last.

A case is defined by the files its reference holds. `cases.manifest()` sorts
them into categories, and each category has one comparison method:

    snapshot  fault.0*.nc        netcdf, per fault, max|diff| must be 0
    restart   fault.r.nc         netcdf
    onfault   fltst_strk*        full series diff
    offfault  srfst_strk*        full series diff
    global    global.dat/.csv    full series diff
    profile   *_strike/_depth    full series diff  (BP8 only writes these)
    cplot     cplot_EQquasi.*    full series diff
    scalars   summary.json       the named values that reference calls important

Nothing here names a benchmark to decide what to check, so a new benchmark is a
row in `cases.CASES` plus a reference directory.

    pytest -m e2e_fast tests/    what CI runs: 101-step cases, BP8's 30 days,
                                 and the clean build
    pytest -m e2e tests/         adds the full BP5 cycle
"""

import json
import os
import subprocess

import pytest

from conftest import ROOT
import cases as C

pytestmark = pytest.mark.e2e

WORK_ROOT = os.path.join(str(ROOT), "work")

ALL = [pytest.param(c, marks=pytest.mark.e2e_fast) if c[4] == "fast"
       else pytest.param(c) for c in C.CASES]
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
        ok = abs(got - want) <= abs(want) * rel
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
