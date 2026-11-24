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
import re
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


@pytest.mark.e2e_fast
def test_mat_het_constant_field_equals_uniform_material():
    """par.mat_het (v1.22.0) feeds a position-dependent shear modulus to the
    element stiffness (src/meshgen.f90), the fault's radiation damping and
    creeping-zone traction (src/faulting.f90) and the CFL step
    (src/read_input.f90). A field that is constant everywhere -- wavelengths
    of 1e15 m, a box larger than the model -- must therefore reproduce a
    uniform model with that same mu bit for bit. It caught a CFL step still
    using the background vp: 1e-3 relative error after 60 steps.
    """
    netCDF4 = pytest.importorskip("netCDF4")
    compset = "test.bp5.qdc.2000"
    tail = {
        "const": "par.nstep = 60\n"
                 "par.mat_het = dict(amp=4e9, lx=1e15, lz=1e15, z0=0.0, "
                 "box=(1e9, 1e9, 1e9))\n",
        "unif": "par.nstep = 60\n"
                "_mu = par.rou*par.vs**2 + 4e9\n"
                "_r = par.vp/par.vs\n"
                "par.vs = (_mu/par.rou)**0.5\n"
                "par.vp = _r*par.vs\n",
    }
    out = {}
    for name, extra in tail.items():
        workdir = os.path.join(WORK_ROOT, f"e2e.mat_het_{name}")
        e = C.env()
        shutil.rmtree(workdir, ignore_errors=True)
        subprocess.run([str(ROOT / "scripts" / "create.newcase"), workdir, compset],
                       cwd=str(ROOT), env=e, check=True, capture_output=True)
        with open(os.path.join(workdir, "user_defined_params.py"), "a") as f:
            f.write(extra)
        subprocess.run(["./case.setup"], cwd=workdir, env=e, check=True,
                       capture_output=True)
        assert os.path.exists(os.path.join(workdir, "input", "mat_het.txt")) == (name == "const")
        r = subprocess.run(["bash", "run.sh"], cwd=workdir, env=e,
                           capture_output=True, text=True)
        assert r.returncode == 0, f"{name}: run.sh exit {r.returncode}\n{r.stdout[-3000:]}"
        out[name] = os.path.join(workdir, "result", "cycle0", "fault.00060.nc")

    with netCDF4.Dataset(out["const"]) as a, netCDF4.Dataset(out["unif"]) as b:
        for v in ("shear_strike", "slip_rate", "state_variable"):
            assert (a[v][:] == b[v][:]).all(), f"{v} differs between a constant mat_het field and the uniform model"


@pytest.mark.e2e_fast
def test_xi_gamma_per_node_step_runs_and_matches_eq15():
    """par.xi_gamma (v1.23.0) runs the per-node time-step factor of Lapusta et
    al. (2000), eq. 15, in src/faulting.f90 -- the real Fortran path, not a
    source-string pin (PR #72 audit, MAJOR).

    test.bp5.qdc.2000 with the cap par.xi raised to 0.5 (its own 0.015 sits
    below every node's eq. 15 value, so the cap would mask the feature) is run
    for 3 steps with xi_gamma = 0 and with xi_gamma = pi/4 from the same state.
    Step 2 is the first step past the base dt: the constant run takes
    dt = 0.5 L/V, the per-node run dt = xi_node L/V at the governing node, so
    their ratio must equal min over velocity-weakening nodes of
    eq. 15 / 0.5 (here 0.18-0.20), up to EQquasi's rounding of dt down to a
    multiple of the base step. A per-node path that is off gives ratio 1.

    Then xi_gamma = 1e-4 makes every velocity-weakening cell larger than h*:
    the run must stop with code 13 and its message, within the timeout, i.e.
    through MPI_Abort rather than a rank-local stop that leaves the others hung.
    """
    import math
    import sys
    import numpy as np
    from defaultParameters import lapusta_xi
    compset, cap, gamma = "test.bp5.qdc.2000", 0.5, math.pi / 4

    def run(name, extra, timeout=900):
        workdir = os.path.join(WORK_ROOT, f"e2e.xi_gamma_{name}")
        e = C.env()
        shutil.rmtree(workdir, ignore_errors=True)
        subprocess.run([str(ROOT / "scripts" / "create.newcase"), workdir, compset],
                       cwd=str(ROOT), env=e, check=True, capture_output=True)
        with open(os.path.join(workdir, "user_defined_params.py"), "a") as f:
            f.write(f"par.nstep = 3\npar.xi = {cap}\n" + extra)
        subprocess.run(["./case.setup"], cwd=workdir, env=e, check=True, capture_output=True)
        r = subprocess.run(["bash", "run.sh"], cwd=workdir, env=e, capture_output=True,
                           text=True, timeout=timeout)
        return workdir, r

    w_off, r_off = run("off", "")
    w_on, r_on = run("on", f"par.xi_gamma = {gamma!r}\n")
    for w, r in ((w_off, r_off), (w_on, r_on)):
        assert r.returncode == 0, f"{w}: run.sh exit {r.returncode}\n{r.stdout[-3000:]}"
    assert "per-node xi" in r_on.stdout and "per-node xi" not in r_off.stdout

    t = {k: np.loadtxt(os.path.join(w, "result", "cycle0", "global.dat"), usecols=(0,))
         for k, w in (("off", w_off), ("on", w_on))}
    dt = {k: np.diff(v) for k, v in t.items()}
    base = min(dt["off"].min(), dt["on"].min())
    assert dt["off"][0] > 2 * base, "step 2 of the constant run is still at the base dt"
    ratio = dt["on"][0] / dt["off"][0]

    # eq. 15 at every velocity-weakening node of this compset (uniform mu, sigma, b).
    sys.path.insert(0, w_on)
    from user_defined_params import par  # the case's own parameters
    sys.path.remove(w_on)
    v = par.on_fault_vars[0]
    mu = par.rou * par.vs ** 2
    xs = [lapusta_xi(gamma * mu / par.dx * L, a * abs(s), b * abs(s), cap)
          for a, b, L, s in zip(v[..., 9].ravel(), v[..., 10].ravel(), v[..., 11].ravel(), v[..., 7].ravel())
          if a < b]
    tol = 2 * base / dt["off"][0]
    assert min(xs) / cap - tol <= ratio <= max(xs) / cap + tol, (
        f"per-node/constant step ratio {ratio:.4f} outside eq. 15 range "
        f"[{min(xs) / cap:.4f}, {max(xs) / cap:.4f}] +- {tol:.4f}")

    w_bad, r_bad = run("abort", "par.xi_gamma = 1e-4\n", timeout=600)
    assert r_bad.returncode != 0, "xi_gamma = 1e-4 (cells larger than h*) did not fail"
    assert "lapusta_xi: no stable step" in r_bad.stdout + r_bad.stderr


def _binary_has_petsc():
    exe = C.installed_binary()
    if exe is None:
        pytest.fail("no built bin/eqquasi-<version>: build it (rule 3a); "
                    "a missing binary is a failure, not a skip")
    r = subprocess.run(["ldd", str(exe)], capture_output=True, text=True)
    return "libpetsc" in r.stdout


# Text outputs. A station file (fltst_, srfst_) is written by the rank that
# owns the station; cplot_, p1output and p2output by every rank, each its own
# fault nodes' lines; the rest by rank 0 from scalars every rank holds.
# compare_series takes each column's tolerance from its printed precision:
# the 7-digit columns (E15.7: srfst_ always; fltst_ and global.dat in the
# bp 8 layout) are compared at about 2e-6, not at DD_RTOL, which applies to
# the columns printed in full double precision (cplot_, p1/p2output,
# peak_sliprate_per_fault, and fltst_ and global.dat otherwise).
DD_TEXT_OUTPUTS = ("fltst_", "srfst_", "cplot_", "p1output", "p2output", "global",
                   "peak_sliprate_per_fault")
# bp 8, the section 4.3 line files (slip_2_, slip_3_, shear_stress_2_,
# shear_stress_3_ and pore_pressure_, each along strike and along depth):
# every rank records the line nodes it owns, rank 0 writes the sum.
DD_BP8_LINES = ("slip_", "shear_stress_", "pore_pressure_")
# The kinds of text output solver 2 writes for each compset of
# test_distributed_loop_matches_petsc_loop, listed from its 31-step run. A
# kind listed here that solver 2 did not write fails the test: a comparison
# over no files passes whatever solver 3 does (rule 2).
DD_EXPECTED_KINDS = {
    "test.bp5.qdc.2000": ("fltst_", "srfst_", "cplot_", "p1output", "p2output", "global", "peak_sliprate_per_fault"),
    "test.bp1002.qdc.zone.2500": ("fltst_", "srfst_", "cplot_", "p1output", "p2output", "global", "peak_sliprate_per_fault"),
    "test.bp7.qdc.a.10": ("fltst_", "srfst_", "cplot_", "p1output", "p2output", "global", "peak_sliprate_per_fault"),
    "test.bp8.qdc.gs.10": ("fltst_", "srfst_", "cplot_", "global", "peak_sliprate_per_fault", "slip_", "shear_stress_", "pore_pressure_"),
}
# The whole-mesh numbering solver 3 computes in closed form, each rank sending
# its own part: integers and grid coordinates, so equal arrays and no tolerance.
DD_MESH_FILES = ("eqquasi.mesh.coor.nc", "eqquasi.mesh.ien.nc", "eqquasi.mesh.nsmp.nc")
# disp.*.nc and fault.*.nc are stored as float64 (the solver's dp) since
# 1.24.8 and are compared at DD_RTOL like the text columns (they were float32,
# compared at 1e-6, before). A node written under another node's number moves
# the field far past it.
# The double-precision columns. Measured 2026-11-13 over 101 steps at 1, 4
# and 8 ranks: worst 1.6e-13 of a column's scale. 1e-10 leaves room for a
# different MUMPS ordering.
DD_RTOL = 1e-10
# Round-off floors of the fltst_ out-of-plane columns (owner decision, board
# row 1, 2026-11-18/19), in the un-logged units of each column: slip_3 in m,
# slip_rate_3 in m/s (the file holds its log10), shear_stress_3 in MPa. An
# entry whose two values differ by no more than the floor passes, and so does
# a log10 entry whose un-logged values both sit under it. The floors apply to
# every fltst_ comparison, the gold comparison included; no srfst_ floor, and
# no other full-cycle comparison rule (srfst_, global.dat, tolerances) changed.
# The constants live in cases.py, where compare_series applies them.
DD_STATION_FLOORS = C.STATION_FLOORS
# The largest rank's share of the whole mesh's nodes at 4 ranks. A block is a
# quarter of the cells plus one ghost line per cut direction; a rank that
# held the whole mesh would be 1.
# Measured 2026-11-13 at 4 ranks: 0.267 and 0.264 of the nodes, 0.25 of the
# elements.
DD_MAX_NODE_SHARE = 0.35


def _assert_distributed_mesh(stdout, nranks):
    """The run summary of par.solver = 3 must show the distributed loop on
    `nranks` ranks, one block each, and no rank holding the whole mesh."""
    def ints(label):
        m = re.search(re.escape(label) + r"\s*=\s*((?:\d+\s+)+)", stdout)
        assert m, f"solver 3 did not print '{label}': not the distributed loop\n{stdout[-2000:]}"
        return [int(v) for v in m.group(1).split()]

    assert ints("= MPI ranks") == [nranks]
    blocks = ints("= Blocks (x, y, z)")
    assert len(blocks) == 3 and blocks[0] * blocks[1] * blocks[2] == nranks, \
        f"blocks {blocks} do not make {nranks} ranks"
    largest, whole = ints("= Nodes on the largest rank, in the mesh")
    assert largest <= DD_MAX_NODE_SHARE * whole, \
        f"a rank holds {largest} of the mesh's {whole} nodes: the mesh is not distributed"
    largest_el, whole_el = ints("= Elements on the largest rank, in the mesh")
    assert largest_el <= DD_MAX_NODE_SHARE * whole_el, \
        f"a rank holds {largest_el} of the mesh's {whole_el} elements"
    return blocks


def _dd_profile_series(gpath, rpath):
    """compare_series for p1output and p2output (time, slip_s, slip_d, tstk,
    tdip; no header, so compare_series scales each column by its own largest
    value). On a fault where the dip traction vanishes by symmetry (bp 7) that
    column is round-off of either solver, about 2e-9 Pa beside 1.6e7 Pa of
    strike traction, and its own maximum is no scale: tdip is judged against
    the vector it is a component of, the larger of its own and tstk's maximum.
    The other columns keep their own scale; every column keeps DD_RTOL."""
    import numpy as np
    g, r = C.read_any(gpath), C.read_any(rpath)
    if g.shape != r.shape:
        return False, f"shape {r.shape} vs reference {g.shape}"
    g, r = np.atleast_2d(np.array(g, float)), np.atleast_2d(np.array(r, float))
    if g.shape[1] != 5:
        return False, f"{g.shape[1]} columns, not time, slip_s, slip_d, tstk, tdip"
    # A non-finite reference entry makes its column's scale, hence its
    # tolerance, infinite: the column would pass whatever the run holds (rule 2).
    if not np.isfinite(g).all():
        i, j = (int(k[0]) for k in np.nonzero(~np.isfinite(g)))
        return False, (f"{os.path.basename(gpath)}: reference column {j+1} holds {g[i, j]} at row {i}; "
                       f"a non-finite reference disables the column's check")
    scale = np.abs(g).max(axis=0)
    scale[4] = max(scale[4], scale[3])
    ok, nbad, idx, worst = C.compare_arrays(g, r, DD_RTOL, DD_RTOL, scale=scale)
    if not ok:
        i, j = idx
        return False, (f"{nbad} of {g.size} entries outside rtol={DD_RTOL:g}; worst at row {i}, "
                       f"column {j+1}: {r[i, j]:.8g} vs {g[i, j]:.8g}")
    return True, f"{g.shape[0]} rows x {g.shape[1]} cols, worst diff = {worst:.1e} of the column scales"


def _dd_station_series(gpath, rpath):
    """compare_series for fltst_ (no header, so compare_series scales each
    column by its own largest value). A dip component that vanishes by
    symmetry (bp 8 at the well node, bp 7) is round-off of either solver and
    its own maximum is no scale: an entry of a dip column also passes when it
    differs by no more than DD_RTOL of the vector it is a component of, the
    largest magnitude among that vector's columns in the same file (slip:
    strike, dip; slip rate: strike or magnitude, dip, compared in m/s where
    the file holds log10; traction: strike, dip, and normal where the file
    holds it). Every other column, and a dip entry above that floor, is
    judged as compare_series judges it: its own scale and its printed
    precision (2e-6 for the 7-digit columns, DD_RTOL otherwise)."""
    import numpy as np
    g, r = C.read_any(gpath), C.read_any(rpath)
    if g.shape != r.shape:
        return False, f"shape {r.shape} vs reference {g.shape}"
    g, r = np.atleast_2d(np.array(g, float)), np.atleast_2d(np.array(r, float))
    # Columns (0-based) by column count, src/library_output.f90: the dip
    # column, the columns of its vector, whether they hold log10.
    layouts = {
        # t, slip_2, slip_3, log10 slip_rate_2, log10 slip_rate_3,
        # shear_stress_2, shear_stress_3, pore pressure, Darcy 2, 3, log10 state
        11: ((2, (1, 2), False, "slip_3"), (4, (3, 4), True, "slip_rate_3"),
             (6, (5, 6), False, "shear_stress_3")),
        # t, slip_s, slip_d, slip rate magnitude, slip rate along dip, tstk,
        # tdip, effective normal stress, log10 state
        9: ((2, (1, 2), False, "slip_3"), (4, (3, 4), False, "slip_rate_3"),
            (6, (5, 6, 7), False, "shear_stress_3")),
    }
    ncols = g.shape[1]
    if ncols not in layouts:
        return False, f"{ncols} columns: not a station layout this comparison knows"
    # A non-finite reference entry makes its column's scale, hence its
    # tolerance, infinite: the column would pass whatever the run holds. No
    # committed reference holds an inf or nan and no case documents one
    # (rg over data/), so there is no both-inf exemption (rule 2).
    if not np.isfinite(g).all():
        i, j = (int(k[0]) for k in np.nonzero(~np.isfinite(g)))
        return False, (f"{os.path.basename(gpath)}: reference column {j+1} holds {g[i, j]} at row {i}; "
                       f"a non-finite reference disables the column's check")
    col_rtol = C.column_printed_rtol(gpath, ncols, DD_RTOL)
    with np.errstate(invalid="ignore", over="ignore"):
        diff = np.abs(r - g)
        own = C._scaled_tol(g, col_rtol, DD_RTOL, np.abs(g).max(axis=0))
        bad = ~(diff <= own)
        for j, vec, log10, name in layouts[ncols]:
            lin = (lambda a: 10.0 ** a) if log10 else (lambda a: a)
            vscale = max(np.abs(lin(g[:, c])).max() for c in vec)
            # ATOL_MACHINE_NOISE (1e-18) is added here, in linear units, to the
            # dip column j only: on the log10 slip-rate column it is 1e-18 m/s,
            # the MPI-reduction noise of a rate that is zero by symmetry (cases.
            # _scaled_tol), nine orders above +-1e-27 and nine below creep.
            # The strike column never enters this `&=`; `own` judges it, in
            # log10 units, where 1e-18 is inert (test_dd_station_series_noise_
            # floor_applies_to_dip_rate_column_only).
            tol = max(DD_RTOL * vscale + C.ATOL_MACHINE_NOISE, DD_STATION_FLOORS[name])
            # No floor on a non-finite raw entry (cases.station_floor_exempt):
            # 10**-inf is 0, which would read a -inf rate as a zero rate.
            finite = np.isfinite(g[:, j]) & np.isfinite(r[:, j])
            bad[:, j] &= ~(finite & (np.abs(lin(r[:, j]) - lin(g[:, j])) <= tol))
        ratio = np.where(bad, diff / np.maximum(own, 1e-300), 0.0)
    if bad.any():
        i, j = np.unravel_index(int(np.argmax(np.nan_to_num(ratio, nan=np.inf))), g.shape)
        cols = sorted({int(c) + 1 for c in np.nonzero(bad)[1]})
        return False, (f"{int(bad.sum())} of {g.size} entries outside rtol={col_rtol[j]:g}, "
                       f"in columns {cols}; worst at row {i}, column {j+1}: "
                       f"{r[i, j]:.8g} vs {g[i, j]:.8g}")
    return True, f"{g.shape[0]} rows x {ncols} cols"


def _dd_differences(s2, s3, expected):
    """Every way solver 3's cycle directory `s3` differs from solver 2's `s2`,
    one line each: the text series, the mesh files and the disp and fault
    snapshots.
    A kind of output in `expected` that solver 2 did not write is itself a
    line (rule 2), and so is a file of a compared kind that only solver 3
    wrote."""
    import numpy as np
    netCDF4 = pytest.importorskip("netCDF4")
    bad = []

    def pair(f):
        if not os.path.exists(os.path.join(s3, f)):
            bad.append(f"{f}: not written by solver 3")
            return None
        return os.path.join(s2, f), os.path.join(s3, f)

    text_kinds = DD_TEXT_OUTPUTS + DD_BP8_LINES
    unknown = [p for p in expected if p not in text_kinds]
    assert not unknown, f"expected kinds {unknown} are not compared kinds {text_kinds}"
    for p in text_kinds:
        files = C.station_files(s2, p)
        if not files and p in expected:
            bad.append(f"{p}*: solver 2 wrote none, nothing to compare")
        if files and p not in expected:
            bad.append(f"{p}*: solver 2 wrote {len(files)}, not a kind expected of this compset")
        for f in files:
            ab = pair(f)
            if ab:
                if p in ("p1output", "p2output"):
                    ok, msg = _dd_profile_series(*ab)
                elif p == "fltst_":
                    ok, msg = _dd_station_series(*ab)
                else:
                    ok, msg = C.compare_series(*ab, rtol=DD_RTOL, floor_rtol=DD_RTOL)
                if not ok:
                    bad.append(f"{f}: {msg}")
        for f in C.station_files(s3, p):
            if f not in files:
                bad.append(f"{f}: written by solver 3 only")

    for f in DD_MESH_FILES:
        if not os.path.exists(os.path.join(s2, f)):
            bad.append(f"{f}: not written by solver 2, nothing to compare")
            continue
        ab = pair(f)
        if not ab:
            continue
        with netCDF4.Dataset(ab[0]) as a, netCDF4.Dataset(ab[1]) as b:
            if set(a.variables) != set(b.variables):
                bad.append(f"{f}: variables {sorted(b.variables)} vs {sorted(a.variables)}")
                continue
            for v in a.variables:
                x, y = np.asarray(a[v][:]), np.asarray(b[v][:])
                if x.shape != y.shape:
                    bad.append(f"{f}: {v} has shape {y.shape} vs {x.shape}")
                elif not np.array_equal(x, y):
                    bad.append(f"{f}: {v} differs in {int((x != y).sum())} of {x.size} entries")

    for kind in ("disp.", "fault."):
        snaps = sorted(f for f in os.listdir(s2) if f.startswith(kind) and f.endswith(".nc"))
        if not snaps:
            bad.append(f"{kind}*.nc: solver 2 wrote none, nothing to compare")
        for f in sorted(os.listdir(s3)):
            if f.startswith(kind) and f.endswith(".nc") and f not in snaps:
                bad.append(f"{f}: written by solver 3 only")
        for f in snaps:
            ab = pair(f)
            if ab:
                ok, msg = C.compare_netcdf(*ab, rtol=DD_RTOL, floor_rtol=DD_RTOL)
                if not ok:
                    bad.append(f"{f}: {msg}")
    return bad


@pytest.mark.e2e_fast
@pytest.mark.parametrize("compset", sorted(DD_EXPECTED_KINDS))
def test_distributed_loop_matches_petsc_loop(compset):
    """par.solver = 3 (each rank owns a block of the mesh, src/
    solveTimeLoopPETScDD.f90) solves the same system with the same MUMPS
    Cholesky as par.solver = 2, so the two must agree to round-off. Run at 4
    ranks -- blocks (2, 2, 1), which cut every fault along strike and put the
    block boundary on the fault plane -- against solver 2 at the compset's own
    rank count, 31 steps each. The two-fault compset is here for rule 12, bp 7
    for a fault whose dip components vanish by symmetry, bp 8 for the pore
    pressure every rank advances on its own fault nodes.
    Compared (_dd_differences): the text series, the mesh files in whole-mesh
    numbering, and the disp and fault snapshots -- every file solver 3 writes
    rank by rank.

    Not in e2e_ci: CI builds -m ubuntu, which has no PETSc.
    """
    if not _binary_has_petsc():
        pytest.skip("this binary was built without PETSc (par.solver = 2, 3 unavailable)")
    cyc, out = _dd_run_pair(compset, compset, "par.solver = 2\npar.nstep = 31\n",
                            "par.solver = 3\npar.nstep = 31\npar.HPC_ncpu = 4\n")
    _assert_distributed_mesh(out, nranks=4)
    bad = _dd_differences(cyc["s2"], cyc["s3"], DD_EXPECTED_KINDS[compset])
    assert not bad, "solver 3 differs from solver 2:\n" + "\n".join(bad)


def _dd_run_one(name, tag, extra):
    """One workflow run of `tag` = (compset, tag) with the parameter lines
    `extra` appended, in scratch/e2e.dd_<name>.<tag>; returns its cycle
    directory and stdout."""
    compset, tag = tag
    workdir = os.path.join(WORK_ROOT, f"e2e.dd_{name}.{tag}")
    e = C.env()
    shutil.rmtree(workdir, ignore_errors=True)
    subprocess.run([str(ROOT / "scripts" / "create.newcase"), workdir, compset],
                   cwd=str(ROOT), env=e, check=True, capture_output=True)
    with open(os.path.join(workdir, "user_defined_params.py"), "a") as f:
        f.write(extra)
    subprocess.run(["./case.setup"], cwd=workdir, env=e, check=True,
                   capture_output=True)
    r = subprocess.run(["bash", "run.sh"], cwd=workdir, env=e,
                       capture_output=True, text=True)
    assert r.returncode == 0, f"{name}: run.sh exit {r.returncode}\n{r.stdout[-3000:]}"
    return os.path.join(workdir, "result", "cycle0"), r.stdout


def _dd_run_pair(compset, tag, s2, other, second="s3"):
    """Run `compset` twice through the workflow, with the parameter lines `s2`
    and `other` appended, in scratch/e2e.dd_s2.<tag> and
    scratch/e2e.dd_<second>.<tag> (second = "s3" for solver 3, "s4" for solver 4).
    Returns the two cycle directories, keyed "s2" and `second`, and the second
    run's stdout."""
    cyc, out = {}, {}
    for name, extra in {"s2": s2, second: other}.items():
        cyc[name], out[name] = _dd_run_one(name, (compset, tag), extra)
    return cyc, out[second]


# Layouts the 4-rank test above does not reach (its blocks are (2, 2, 1)):
# id -> compset, ranks of the solver 3 run, parameter lines for both runs.
# MPI_Dims_create gives (2, 2, 2) at 8 ranks and (3, 2, 2) at 12, the 3 in
# the direction with the most cells. Solver 2 runs at 4 ranks: its result
# does not depend on the rank count, and the host is shared (rule 15).
DD_LAYOUT_CASES = {
    "bp1002.12": ("test.bp1002.qdc.zone.2500", 12, ""),
    "bp7.8": ("test.bp7.qdc.a.10", 8, ""),
    "bp8.gs.8": ("test.bp8.qdc.gs.10", 8, ""),
    # The Peaceman well of BP8-QD-PW on the GS test mesh, with Table 1's Q0
    # for PW. The well node (0, 0) is the grid point where the four in-plane
    # blocks meet.
    "bp8.pw.8": ("test.bp8.qdc.gs.10", 8, "par.fluid_src = 2\npar.fluid_q0 = 0.0015\n"),
}


@pytest.mark.e2e_fast
@pytest.mark.parametrize("case", sorted(DD_LAYOUT_CASES))
def test_distributed_loop_matches_petsc_loop_z_and_corner_blocks(case):
    """test_distributed_loop_matches_petsc_loop on the layouts it leaves out
    (PR #83 audit): more than one block along z, so the halo exchange crosses
    z; more than one block along x and along z, so four blocks meet at a line
    through the fault plane (y is its normal) and a fault node there has
    diagonal neighbours on another rank; bp 8 with its section 4.3 depth line
    (z from -400 to 400 m at x = 0), which a z boundary at z = 0 would split
    between ranks (expected from the geometry only: the code does not print
    which ranks record its slots and nothing here asserts it); and the
    fluid_src = 2 well.
    Same comparison and tolerances as the 4-rank test. The layout is read
    from the run summary: a run that fell back to one block in z fails.
    """
    if not _binary_has_petsc():
        pytest.skip("this binary was built without PETSc (par.solver = 2, 3 unavailable)")
    compset, nranks, extra = DD_LAYOUT_CASES[case]
    cyc, out = _dd_run_pair(compset, case,
                            "par.solver = 2\npar.nstep = 31\npar.HPC_ncpu = 4\n" + extra,
                            f"par.solver = 3\npar.nstep = 31\npar.HPC_ncpu = {nranks}\n" + extra)
    blocks = _assert_distributed_mesh(out, nranks=nranks)
    print(f"\n{case}: {nranks} ranks, blocks (x, y, z) = {tuple(blocks)}")
    assert blocks[2] > 1, f"blocks {blocks}: one block along z, no z halo exchange"
    assert blocks[0] > 1, f"blocks {blocks}: one block along x, no in-plane block corner"
    if "fluid_src = 2" in extra:
        assert "well index" in out, "fluid_src = 2 did not report its well index: not the well path"
    bad = _dd_differences(cyc["s2"], cyc["s3"], DD_EXPECTED_KINDS[compset])
    assert not bad, "solver 3 differs from solver 2:\n" + "\n".join(bad)


DD_KILL_GRACE = 30   # s a SIGKILLed process group has to close its stdout


def _kill_group_and_collect(p, grace):
    """SIGKILL p's process group and return what it printed. A member that left
    the group (its own session) keeps the pipe open and communicate() would
    block forever: past `grace` seconds that is a failure naming the pid."""
    import signal
    os.killpg(p.pid, signal.SIGKILL)
    try:
        out, _ = p.communicate(timeout=grace)
    except subprocess.TimeoutExpired:
        pytest.fail(f"process group {p.pid} still holds the output pipe {grace} s after SIGKILL: "
                    f"a member left the group and is still running")
    return out


@pytest.mark.e2e_fast
def test_distributed_abort_on_one_rank_stops_every_rank():
    """check_effective_normal (src/faulting.f90) is reached only by the rank
    that owns the unclamped fault node. Under par.solver = 3 a bare `stop`
    there would leave the other ranks blocked in the next collective; it
    calls MPI_ABORT instead. test.bp8.qdc.gs.10 at 8 ranks (the run summary
    with the block layout is printed at the end of a run, so an aborted run
    has none and the layout is not asserted here), with Q0 raised from 0.003 to 1 m^3/s: the pore pressure at the injection
    node passes the 25 MPa of initial effective normal stress in the first
    steps. The whole mpirun must then exit nonzero within the timeout, with
    the routine's own message and MPI_ABORT's error code 508 in the output.
    A run that ends cleanly, ends with another message or code, or hangs
    fails.
    """
    if not _binary_has_petsc():
        pytest.skip("this binary was built without PETSc (par.solver = 3 unavailable)")
    compset, nranks = "test.bp8.qdc.gs.10", 8
    workdir = os.path.join(WORK_ROOT, "e2e.dd_abort")
    e = C.env()
    shutil.rmtree(workdir, ignore_errors=True)
    subprocess.run([str(ROOT / "scripts" / "create.newcase"), workdir, compset],
                   cwd=str(ROOT), env=e, check=True, capture_output=True)
    with open(os.path.join(workdir, "user_defined_params.py"), "a") as f:
        f.write(f"par.solver = 3\npar.HPC_ncpu = {nranks}\npar.fluid_q0 = 1.0\n")
    subprocess.run(["./case.setup"], cwd=workdir, env=e, check=True, capture_output=True)
    # Its own process group: on a hang, every rank is killed with the launcher.
    p = subprocess.Popen(["bash", "run.sh"], cwd=workdir, env=e, stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT, text=True, start_new_session=True)
    try:
        out, _ = p.communicate(timeout=600)
    except subprocess.TimeoutExpired:
        out = _kill_group_and_collect(p, DD_KILL_GRACE)
        pytest.fail(f"par.solver = 3 did not stop within 600 s: a rank-local abort hangs the job\n{out[-3000:]}")
    # run.sh echoes only the tail of a failed cycle's log; the summary and
    # the message are in the log itself.
    with open(os.path.join(workdir, "log", "cycle0.log")) as f:
        out += f.read()
    banners = out.count("EFFECTIVE NORMAL STRESS IS NO LONGER COMPRESSIVE")
    print(f"\nabort: exit {p.returncode}, the message appears {banners} time(s) in run.sh output plus log")
    assert banners >= 1, f"check_effective_normal was not reached (exit {p.returncode}):\n{out[-3000:]}"
    assert p.returncode != 0, "the run printed the abort message and still exited 0"
    assert "MPI_ABORT" in out and re.search(r"[Ee]rrorcode:?\s*508\b", out), \
        f"no MPI_ABORT with error code 508 in the output (exit {p.returncode}):\n{out[-3000:]}"


# --- par.solver = 4, slip to traction: each solver-3 case above, forced to 4 ---
# Steps 2..30 of the 31-step runs are the operator's matrix-vector product alone;
# step 1 (a disp snapshot) and step 31 (the stop step) also solve in full.
S2T_NSTEP = 31
# The operator against the solve inside the run (src/solveTimeLoopPETScDD.f90,
# ddS2tStations prints it on every full-solve step). The bound is the solver's own:
# 1e-8 of each field's largest value, floored at 1e-14 of the operator's largest
# entry (maxabs_S for the traction, maxabs_Su for the displacement, both printed
# on the check line). The floor is the owner's (2026-11-19, "go").
S2T_CHECK_RTOL = 1e-8
S2T_CHECK_FLOOR = 1e-14
S2T_CHECK_RE = (r"EQ_S2T check step maxdiff_traction maxabs_traction maxdiff_u maxabs_u maxabs_S maxabs_Su\s+"
                r"(\d+)" + r"\s+(\S+)" * 6)


def _s2t_violations(rows):
    """The check rows (step, then six floats) that exceed the solver's bound.
    This is a Python copy of the bound the solver applies (src/solveTimeLoopPETScDD.f90,
    max(S2T_CHECK_RTOL*field_max, S2T_CHECK_FLOOR*operator_max)). The Fortran
    abort on a wrong operator is tested by test_slip_to_traction_planted_operator_row_aborts."""
    bad = []
    for r in rows:
        dt_, at_, du, au, st, su = (float(v) for v in r[1:])
        bt = max(S2T_CHECK_RTOL * at_, S2T_CHECK_FLOOR * st)
        bu = max(S2T_CHECK_RTOL * au, S2T_CHECK_FLOOR * su)
        if not (dt_ <= bt and du <= bu):
            bad.append(f"step {r[0]}: operator vs solve, traction {dt_:g} (bound {bt:g}), "
                       f"displacement {du:g} (bound {bu:g})")
    return bad


def _assert_s2t_ran(cyc, nstep):  # nstep: the run's last step
    """The run of `cyc` (a cycle directory) built the operator, and its
    full-solve steps hold the first and the last step, each with the operator's
    traction and station displacement within S2T_CHECK_RTOL of the solve's.
    A run that stayed on the solve (no EQ_S2T lines) fails (rule 2)."""
    log_path = os.path.join(os.path.dirname(os.path.dirname(cyc)), "log", "cycle0.log")
    with open(log_path) as f:
        log = f.read()
    assert re.search(r"EQ_S2T built sec ncol nbat\s+\S+\s+\d+\s+\d+", log), \
        f"solver 4 did not print its operator build line: the slip-to-traction path did not run\n{log[-2000:]}"
    rows = re.findall(S2T_CHECK_RE, log)
    steps = [int(r[0]) for r in rows]
    assert 1 in steps and nstep in steps and len(steps) < nstep, \
        f"full-solve steps {steps}: want step 1 and step {nstep} and fewer than {nstep} of them (else the operator is never applied alone)"
    bad = _s2t_violations(rows)
    assert not bad, "\n".join(bad)


@pytest.mark.e2e_fast
def test_s2t_check_assert_catches_a_planted_wrong_operator():
    """Checks the Python copy of the bound (_s2t_violations), not the Fortran
    abort (that is test_slip_to_traction_planted_operator_row_aborts).
    The assert of _assert_s2t_ran on planted check rows, no run.
    Largest traction 1e6, largest S entry 1e3: relative bound 1e-8 * 1e6 = 1e-2,
    floor 1e-14 * 1e3 = 1e-11. A traction difference of 1e-5 relative (10) must
    fail, 1e-9 relative (1e-3, below the bound) must pass; with a traction zero
    by physics (largest 0) a difference of 1e-12 (below the floor) passes and 1e-9
    (above it) fails; a displacement difference of 1e-5 against largest 1 fails."""
    f = lambda *v: tuple(f"{x:e}" for x in v)
    assert _s2t_violations([("1",) + f(1e6 * 1e-5, 1e6, 0.0, 1.0, 1e3, 1.0)]), "1e-5 relative not caught"
    assert not _s2t_violations([("1",) + f(1e6 * 1e-9, 1e6, 0.0, 1.0, 1e3, 1.0)]), "1e-9 relative rejected"
    assert not _s2t_violations([("1",) + f(1e-12, 0.0, 0.0, 1.0, 1e3, 1.0)]), "floor not applied"
    assert _s2t_violations([("1",) + f(1e-9, 0.0, 0.0, 1.0, 1e3, 1.0)]), "above the floor not caught"
    assert _s2t_violations([("1",) + f(0.0, 1.0, 1.0e-5, 1.0, 1e3, 1.0)]), "displacement 1e-5 not caught"


# The solver's own abort on a wrong operator (src/solveTimeLoopPETScDD.f90,
# ddS2tStations), and its refusal of a row outside the operator.
S2T_PLANT_ABORT = ("EQ_S2T: operator disagrees with the solve beyond max(1e-8 of the field maximum, "
                   "1e-14 of the operator maximum)")
S2T_PLANTS = {
    # row 0 doubled after the build: step 1 solves in full and must stop the run.
    "row0": ("-eq_s2t_plant_row", 0, S2T_PLANT_ABORT),
    "out_of_range": ("-eq_s2t_plant_row", 1000000000,
                     "EQ_S2T: -eq_s2t_plant_row 1000000000 is outside the operator rows [0, "),
    # row 0 set to NaN after the build: the finite check of the operator product
    # S g of step 1 (which runs before the operator check's NaN flag) must stop the run.
    "nan_row0": ("-eq_s2t_plant_nan", 0, "EQ_S2T: NaN or Inf in the operator product S g, step"),
    "nan_out_of_range": ("-eq_s2t_plant_nan", 1000000000,
                         "EQ_S2T: -eq_s2t_plant_nan 1000000000 is outside the operator rows [0, "),
}


@pytest.mark.e2e_fast
@pytest.mark.parametrize("case", sorted(S2T_PLANTS))
def test_slip_to_traction_planted_operator_row_aborts(case):
    """par.solver = 4 with PETSC_OPTIONS -eq_s2t_plant_row <row>: the solver
    doubles that global row of the operator after the build, and the operator
    check of the first full-solve step (step 1) stops the run with a nonzero
    exit and the check's abort message; a row outside [0, n) stops it with the
    hook's own message before any step. -eq_s2t_plant_nan <row> sets that row to
    NaN instead: the finite check of the operator product S g stops the run at step 1.
    31 steps, 4 ranks, test.bp5.qdc.2000."""
    if not _binary_has_petsc():
        pytest.skip("this binary was built without PETSc (par.solver = 4 unavailable)")
    option, row, message = S2T_PLANTS[case]
    workdir = os.path.join(WORK_ROOT, "e2e.dd_plant_s4." + case)
    e = C.env()
    e["PETSC_OPTIONS"] = f"{option} {row}"
    shutil.rmtree(workdir, ignore_errors=True)
    subprocess.run([str(ROOT / "scripts" / "create.newcase"), workdir, "test.bp5.qdc.2000"],
                   cwd=str(ROOT), env=e, check=True, capture_output=True)
    with open(os.path.join(workdir, "user_defined_params.py"), "a") as f:
        f.write(f"par.solver = 4\npar.nstep = {S2T_NSTEP}\npar.HPC_ncpu = 4\n")
    subprocess.run(["./case.setup"], cwd=workdir, env=e, check=True, capture_output=True)
    r = subprocess.run(["bash", "run.sh"], cwd=workdir, env=e, capture_output=True, text=True, timeout=900)
    log_path = os.path.join(workdir, "log", "cycle0.log")
    out = r.stdout + r.stderr + (open(log_path).read() if os.path.exists(log_path) else "")
    assert r.returncode != 0, f"{case}: run.sh exit 0 with a planted operator row\n{out[-3000:]}"
    assert "EQ_S2T built sec" in out, f"{case}: the operator was not built\n{out[-3000:]}"
    assert message in out, f"{case}: exit {r.returncode} without '{message}'\n{out[-3000:]}"
    if case == "row0":
        assert "EQ_S2T: -eq_s2t_plant_row 0: that operator row is doubled" in out, out[-3000:]
    if case == "nan_row0":
        assert "EQ_S2T: -eq_s2t_plant_nan 0: that operator row is set to NaN" in out, out[-3000:]
        # The abort prints `it` list-directed (gfortran pads the integer), so
        # match "step" then whitespace then exactly 1.
        assert re.search(r"operator product S g, step\s+1\s*$", out, re.M), \
            f"the NaN abort did not name step 1\n{out[-3000:]}"


@pytest.mark.e2e_fast
@pytest.mark.parametrize("compset", sorted(DD_EXPECTED_KINDS))
def test_slip_to_traction_matches_petsc_loop(compset):
    """test_distributed_loop_matches_petsc_loop with par.solver = 4: same
    compsets, 4 ranks, 31 steps, same comparison at the same DD_RTOL against
    solver 2, and the operator's own check lines (_assert_s2t_ran)."""
    if not _binary_has_petsc():
        pytest.skip("this binary was built without PETSc (par.solver = 2, 4 unavailable)")
    cyc, out = _dd_run_pair(compset, "s4." + compset, "par.solver = 2\npar.nstep = 31\n",
                            "par.solver = 4\npar.nstep = 31\npar.HPC_ncpu = 4\n", second="s4")
    _assert_distributed_mesh(out, nranks=4)
    _assert_s2t_ran(cyc["s4"], S2T_NSTEP)
    bad = _dd_differences(cyc["s2"], cyc["s4"], DD_EXPECTED_KINDS[compset])
    assert not bad, "solver 4 differs from solver 2:\n" + "\n".join(bad)


@pytest.mark.e2e_fast
@pytest.mark.parametrize("case", sorted(DD_LAYOUT_CASES))
def test_slip_to_traction_matches_petsc_loop_z_and_corner_blocks(case):
    """test_distributed_loop_matches_petsc_loop_z_and_corner_blocks with
    par.solver = 4: the same layouts (more than one block along z and along x,
    the bp 8 well), the same comparison."""
    if not _binary_has_petsc():
        pytest.skip("this binary was built without PETSc (par.solver = 2, 4 unavailable)")
    compset, nranks, extra = DD_LAYOUT_CASES[case]
    cyc, out = _dd_run_pair(compset, "s4." + case,
                            "par.solver = 2\npar.nstep = 31\npar.HPC_ncpu = 4\n" + extra,
                            f"par.solver = 4\npar.nstep = 31\npar.HPC_ncpu = {nranks}\n" + extra,
                            second="s4")
    blocks = _assert_distributed_mesh(out, nranks=nranks)
    assert blocks[2] > 1 and blocks[0] > 1, f"blocks {blocks}: no z halo exchange or no in-plane block corner"
    _assert_s2t_ran(cyc["s4"], S2T_NSTEP)
    bad = _dd_differences(cyc["s2"], cyc["s4"], DD_EXPECTED_KINDS[compset])
    assert not bad, "solver 4 differs from solver 2:\n" + "\n".join(bad)


@pytest.mark.e2e_fast
def test_slip_to_traction_stop_on_a_step_that_is_not_an_output_step():
    """A run that ends on par.fluid_tend (bp 8) stops on a step that is neither
    the last of par.nstep nor one of the nt_out output steps, so only the
    `stoptag == 1` term of the full-solve condition (src/solveTimeLoopPETScDD.f90,
    s2tFull) makes the stop step solve in full. Without it the stop step's disp
    snapshot would hold the previous solve's displacement.

    The stop time is the midpoint of two step times of a solver 2 probe run;
    solver 2 and solver 4 then stop on the same step, the solver 4 log holds a
    full-solve line for it, and every file, the disp snapshot of the stop step
    included, agrees at DD_RTOL (_dd_differences)."""
    if not _binary_has_petsc():
        pytest.skip("this binary was built without PETSc (par.solver = 2, 4 unavailable)")
    compset = "test.bp8.qdc.gs.10"
    probe, _ = _dd_run_one("s2", (compset, "s4stop.probe"), "par.solver = 2\npar.nstep = 31\n")
    times = C.read_any(os.path.join(probe, "global.dat"))[:, 0]
    tend = 0.5 * (times[13] + times[14])
    extra = f"par.nstep = 31\npar.fluid_tend = {tend!r}\n"
    cyc, out = _dd_run_pair(compset, "s4stop", "par.solver = 2\n" + extra,
                            "par.solver = 4\npar.HPC_ncpu = 4\n" + extra, second="s4")
    snaps = sorted(int(m.group(1)) for n_ in os.listdir(cyc["s2"])
                   if (m := re.fullmatch(r"disp\.(\d{5})\.nc", n_)))
    stop = snaps[-1]
    nt_out = 100
    assert 1 < stop < 31 and stop % nt_out != 1, \
        f"solver 2 stopped on step {stop}: want one that is neither the last (31) nor an output step"
    assert os.path.exists(os.path.join(cyc["s4"], f"disp.{stop:05d}.nc")), \
        f"solver 4 wrote no disp snapshot for the stop step {stop}"
    _assert_s2t_ran(cyc["s4"], stop)
    bad = _dd_differences(cyc["s2"], cyc["s4"], DD_EXPECTED_KINDS[compset])
    assert not bad, "solver 4 differs from solver 2 at the stop step:\n" + "\n".join(bad)


@pytest.mark.e2e_fast
def test_slip_to_traction_abort_on_one_rank_stops_every_rank():
    """test_distributed_abort_on_one_rank_stops_every_rank with par.solver = 4:
    the unclamped fault node's rank aborts inside the apply loop (the operator
    is built first), and every rank must stop with error code 508."""
    if not _binary_has_petsc():
        pytest.skip("this binary was built without PETSc (par.solver = 4 unavailable)")
    compset, nranks = "test.bp8.qdc.gs.10", 8
    workdir = os.path.join(WORK_ROOT, "e2e.dd_abort_s4")
    e = C.env()
    shutil.rmtree(workdir, ignore_errors=True)
    subprocess.run([str(ROOT / "scripts" / "create.newcase"), workdir, compset],
                   cwd=str(ROOT), env=e, check=True, capture_output=True)
    with open(os.path.join(workdir, "user_defined_params.py"), "a") as f:
        f.write(f"par.solver = 4\npar.HPC_ncpu = {nranks}\npar.fluid_q0 = 1.0\n")
    subprocess.run(["./case.setup"], cwd=workdir, env=e, check=True, capture_output=True)
    p = subprocess.Popen(["bash", "run.sh"], cwd=workdir, env=e, stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT, text=True, start_new_session=True)
    try:
        out, _ = p.communicate(timeout=600)
    except subprocess.TimeoutExpired:
        out = _kill_group_and_collect(p, DD_KILL_GRACE)
        pytest.fail(f"par.solver = 4 did not stop within 600 s: a rank-local abort hangs the job\n{out[-3000:]}")
    with open(os.path.join(workdir, "log", "cycle0.log")) as f:
        out += f.read()
    assert "EQ_S2T built sec" in out, f"the operator was not built before the abort:\n{out[-3000:]}"
    assert "EFFECTIVE NORMAL STRESS IS NO LONGER COMPRESSIVE" in out, \
        f"check_effective_normal was not reached (exit {p.returncode}):\n{out[-3000:]}"
    assert p.returncode != 0, "the run printed the abort message and still exited 0"
    assert "MPI_ABORT" in out and re.search(r"[Ee]rrorcode:?\s*508\b", out), \
        f"no MPI_ABORT with error code 508 in the output (exit {p.returncode}):\n{out[-3000:]}"


# Unsupported combinations stop the run with a message; none falls back to another solver.
S2T_REFUSALS = {
    # (parameter lines, PETSC_OPTIONS, message, needs PETSc). The solver 5 guard
    # (src/eqquasi.f90) runs before any PETSc branch, so it is tested on every build.
    "solver5": ("par.solver = 5\n", None, "par.solver must be 1, 2, 3 or 4; got", False),
    # The message names the run's own solver (F5): a wrong value fails both.
    "krylov4": ("par.solver = 4\n", "-ksp_type cg -pc_type jacobi",
                "par.solver=4: only KSPPREONLY + PCCHOLESKY", True),
    "krylov3": ("par.solver = 3\n", "-ksp_type cg -pc_type jacobi",
                "par.solver=3: only KSPPREONLY + PCCHOLESKY", True),
    # A direct factor from another package: the operator is built by
    # MatMatTransposeSolve on the MUMPS factor, so the run must stop, not go on.
    # nonmumps4: PETSc itself refuses at KSPSetUp, before the solver's own check.
    # superlu4 (Cholesky kept, factor package replaced) is accepted by KSPSetUp and
    # reaches the check in ddS2tBuild, which names the type found. With
    # -pc_type lu the earlier PC-type guard stops the run instead.
    "nonmumps4": ("par.solver = 4\n", "-pc_factor_mat_solver_type petsc",
                  "MatSolverType petsc does not support matrix type mpiaij", True),
    "superlu4": ("par.solver = 4\n", "-pc_factor_mat_solver_type superlu_dist",
                 "EQ_S2T: the factor solver type is 'superlu_dist', not mumps", True),
}


@pytest.mark.e2e_fast
@pytest.mark.parametrize("case", sorted(S2T_REFUSALS))
def test_slip_to_traction_unsupported_combination_fails_loudly(case):
    """par.solver = 5, and par.solver = 4 with PETSc options that replace the
    MUMPS Cholesky factor the operator is built on, end with a nonzero exit and
    the refusal's message, and write no cycle result."""
    extra, petsc_options, message, needs_petsc = S2T_REFUSALS[case]
    if needs_petsc and not _binary_has_petsc():
        pytest.skip("this binary was built without PETSc (par.solver = 3, 4 unavailable)")
    workdir = os.path.join(WORK_ROOT, "e2e.dd_refuse." + case)
    e = C.env()
    if petsc_options:
        e["PETSC_OPTIONS"] = petsc_options
    shutil.rmtree(workdir, ignore_errors=True)
    subprocess.run([str(ROOT / "scripts" / "create.newcase"), workdir, "test.bp5.qdc.2000"],
                   cwd=str(ROOT), env=e, check=True, capture_output=True)
    with open(os.path.join(workdir, "user_defined_params.py"), "a") as f:
        f.write(extra + "par.HPC_ncpu = 2\npar.nstep = 3\n")
    subprocess.run(["./case.setup"], cwd=workdir, env=e, check=True, capture_output=True)
    r = subprocess.run(["bash", "run.sh"], cwd=workdir, env=e, capture_output=True, text=True, timeout=600)
    log_path = os.path.join(workdir, "log", "cycle0.log")
    log = open(log_path).read() if os.path.exists(log_path) else ""
    assert r.returncode != 0, f"{case}: run.sh exit 0, the combination was accepted"
    if message is not None:
        assert message in r.stdout + r.stderr + log, \
            f"{case}: exit {r.returncode} without '{message}'\n{(r.stdout + log)[-2000:]}"
    assert "EQ_S2T built sec" not in log, f"{case}: the operator was built before the refusal"


@pytest.mark.e2e_fast
def test_slip_to_traction_refused_on_a_build_without_petsc():
    """par.solver = 4 on a binary built without PETSc ends with a nonzero exit
    and the refusal message (src/eqquasi.f90), never another solver's result.
    The binary under test is the installed one when it lacks PETSc; otherwise a
    copy of src/ is built here with the PETSc define and objects removed, which
    is also the rule 17a link check (only conda-linux builds PETSc, so the
    machine is known)."""
    exe = C.installed_binary()
    if _binary_has_petsc():
        build = os.path.join(WORK_ROOT, "e2e.nopetsc_build")
        shutil.rmtree(build, ignore_errors=True)
        shutil.copytree(ROOT / "src", build)
        e = C.env()
        e["MACHINE"] = "conda-linux"
        e["EQQUASIROOT"] = str(ROOT)
        r = subprocess.run(["make", "PETSC_DEFINE=", "OBJ_PETSC=", "PETSC_LIB="], cwd=build, env=e,
                           capture_output=True, text=True, timeout=1800)
        assert r.returncode == 0, "the build without PETSc failed to link:\n" + (r.stdout + r.stderr)[-3000:]
        exe = os.path.join(build, "eqquasi")
        assert "libpetsc" not in subprocess.run(["ldd", exe], capture_output=True, text=True).stdout
    workdir = os.path.join(WORK_ROOT, "e2e.dd_refuse.nopetsc4")
    e = C.env()
    shutil.rmtree(workdir, ignore_errors=True)
    subprocess.run([str(ROOT / "scripts" / "create.newcase"), workdir, "test.bp5.qdc.2000"],
                   cwd=str(ROOT), env=e, check=True, capture_output=True)
    with open(os.path.join(workdir, "user_defined_params.py"), "a") as f:
        f.write("par.solver = 4\npar.HPC_ncpu = 2\npar.nstep = 3\n")
    subprocess.run(["./case.setup"], cwd=workdir, env=e, check=True, capture_output=True)
    run_sh = os.path.join(workdir, "run.sh")
    text = open(run_sh).read()
    patched = re.sub(r"eqquasi-[0-9][0-9.]* 2>&1", str(exe) + " 2>&1", text)
    assert patched != text, "run.sh names no eqquasi-<version> binary to replace"
    open(run_sh, "w").write(patched)
    r = subprocess.run(["bash", "run.sh"], cwd=workdir, env=e, capture_output=True, text=True, timeout=600)
    log_path = os.path.join(workdir, "log", "cycle0.log")
    log = open(log_path).read() if os.path.exists(log_path) else ""
    assert r.returncode != 0, "run.sh exit 0: par.solver = 4 ran on a build without PETSc"
    assert "needs MACHINE=conda-linux" in r.stdout + r.stderr + log, \
        f"exit {r.returncode} without the refusal message\n{(r.stdout + log)[-2000:]}"
    assert "EQ_S2T built sec" not in log


# --- the comparison's own contracts: milliseconds, no run, part of the sweep ---

def _write_station_file(path, rows):
    """An 11-column fltst_ file (t, slip_2, slip_3, log10 rate_2, log10 rate_3,
    tau_2, tau_3, pore pressure, Darcy 2, 3, log10 state) in the solver's
    seven-significant-figure exponential format."""
    with open(path, "w") as f:
        for row in rows:
            f.write("".join(f"{v:21.13E}" if k == 0 else f"{v:15.6E}" for k, v in enumerate(row)) + "\n")


def _station_row(t, log_rate_strike, log_rate_dip):
    return [t, 1.0, 1.0, log_rate_strike, log_rate_dip, 1.0e6, 1.0e6, 1.0, 1.0, 1.0, -5.0]


def _station_pair(tmp_path, ref_rows, run_rows):
    g, r = str(tmp_path / "ref" / "fltst_a.txt"), str(tmp_path / "run" / "fltst_a.txt")
    for p, rows in ((g, ref_rows), (r, run_rows)):
        os.makedirs(os.path.dirname(p))
        _write_station_file(p, rows)
    return g, r


@pytest.mark.e2e_fast
@pytest.mark.parametrize("ref_value", [float("inf"), float("-inf"), float("nan")])
def test_dd_station_series_non_finite_reference_is_a_failure_naming_file_and_column(tmp_path, ref_value):
    """An inf in a reference column made its scale inf, hence its tolerance
    inf, and the column passed whatever the run held (rule 2). Here the run
    holds a finite, wildly different rate, and the reference is non-finite in
    the dip log10 slip-rate column (column 5)."""
    ref = [_station_row(1.0, -9.0, -9.0), _station_row(2.0, -9.0, ref_value)]
    run = [_station_row(1.0, -9.0, -9.0), _station_row(2.0, -9.0, -3.0)]
    g, r = _station_pair(tmp_path, ref, run)
    ok, msg = _dd_station_series(g, r)
    assert not ok and "fltst_a.txt" in msg and "column 5" in msg, f"ok={ok}, msg={msg!r}"


@pytest.mark.e2e_fast
def test_dd_station_series_equal_finite_files_pass(tmp_path):
    """Control of the tests around it: identical finite files pass."""
    rows = [_station_row(1.0, -9.0, -9.0), _station_row(2.0, -8.0, -9.5)]
    g, r = _station_pair(tmp_path, rows, rows)
    ok, msg = _dd_station_series(g, r)
    assert ok, msg


@pytest.mark.e2e_fast
def test_dd_station_series_noise_floor_applies_to_dip_rate_column_only(tmp_path):
    """Both runs hold a rate of 1e-20 m/s (log10 -20 against -19.9: 2.6e-21 m/s
    apart, under the 1e-12 m/s slip_rate_3 floor). In the dip column (5) that
    is noise and passes; the same pair in the strike column (4) is a 0.1 shift
    of a log10 value, the floor does not reach it, and it fails."""
    ref = [_station_row(1.0, -20.0, -20.0), _station_row(2.0, -20.0, -20.0)]
    dip = [_station_row(1.0, -20.0, -20.0), _station_row(2.0, -20.0, -19.9)]
    strike = [_station_row(1.0, -20.0, -20.0), _station_row(2.0, -19.9, -20.0)]
    ok_dip, msg_dip = _dd_station_series(*_station_pair(tmp_path / "d", ref, dip))
    ok_strike, msg_strike = _dd_station_series(*_station_pair(tmp_path / "s", ref, strike))
    assert ok_dip, f"dip column under the floor must pass: {msg_dip}"
    assert not ok_strike and "columns [4]" in msg_strike, f"strike column must fail: {msg_strike}"


@pytest.mark.e2e_fast
def test_dd_station_series_dip_rate_above_noise_floor_fails(tmp_path):
    """The floor is 1e-12 m/s and no more: -10 against -9.9 is 2.6e-11 m/s."""
    ref = [_station_row(1.0, -20.0, -10.0), _station_row(2.0, -20.0, -10.0)]
    run = [_station_row(1.0, -20.0, -10.0), _station_row(2.0, -20.0, -9.9)]
    ok, msg = _dd_station_series(*_station_pair(tmp_path, ref, run))
    assert not ok and "columns [5]" in msg, f"ok={ok}, msg={msg!r}"


def _floor_case(tmp_path, col, ref_v, run_v, ncols=11):
    """Reference and run rows equal except column `col` (0-based) at row 1.
    The strike column of the same vector is zeroed (log10 rate: left at -9)
    so the vector's scale does not dwarf the floor under test."""
    def row(t, v):
        r = _station_row(t, -9.0, -9.0)
        if col in (2, 6):
            r[col - 1] = 0.0
        r[col] = v
        return r
    ref = [row(1.0, ref_v), row(2.0, ref_v)]
    run = [row(1.0, ref_v), row(2.0, run_v)]
    return _dd_station_series(*_station_pair(tmp_path, ref, run))


# (column, name, reference value, run value under the floor, run value above it)
_FLOOR_CASES = [
    (2, "slip_3", 0.0, 5e-16, 3e-15),
    (4, "slip_rate_3", -20.0, -13.0, -11.0),        # log10: 1e-13 and 1e-11 m/s
    (6, "shear_stress_3", 0.0, 5e-10, 3e-9),
]


@pytest.mark.e2e_fast
@pytest.mark.parametrize("col,name,ref_v,below,above", _FLOOR_CASES)
def test_dd_station_series_difference_under_floor_passes(tmp_path, col, name, ref_v, below, above):
    ok, msg = _floor_case(tmp_path, col, ref_v, below)
    assert ok, f"{name}: a difference under its floor must pass: {msg}"


@pytest.mark.e2e_fast
@pytest.mark.parametrize("col,name,ref_v,below,above", _FLOOR_CASES)
def test_dd_station_series_difference_above_floor_fails(tmp_path, col, name, ref_v, below, above):
    ok, msg = _floor_case(tmp_path, col, ref_v, above)
    assert not ok and f"column {col+1}" in msg, f"{name}: ok={ok}, msg={msg!r}"


@pytest.mark.e2e_fast
def test_dd_station_series_floors_do_not_reach_in_plane_columns(tmp_path):
    """Shear_stress_2 (column 6) is in-plane. Reference 0 and run 5e-10 with
    its sibling shear_stress_3 also 0 is a difference under the 1e-9
    shear_stress_3 floor and above every other tolerance (scale 0): it must
    fail. A floor wrongly applied to column 6 would pass it."""
    ok, msg = _in_plane_below_floor_case(tmp_path)
    assert not ok and "column 6" in msg, f"ok={ok}, msg={msg!r}"


@pytest.mark.e2e_fast
@pytest.mark.parametrize("ref_v,run_v", [(-20.0, float("-inf")), (float("-inf"), -20.0)])
def test_dd_station_series_non_finite_log10_slip_rate_is_never_floored(tmp_path, ref_v, run_v):
    """10**-inf is 0: a -inf log10 rate in the run (or reference) against a
    -20 one is 1e-20 m/s apart in linear units, under the 1e-12 floor, and
    must still fail (cases.station_floor_exempt: no floor on a non-finite
    raw entry)."""
    ok, msg = _floor_case(tmp_path, 4, ref_v, run_v)
    assert not ok and "column 5" in msg, f"ref={ref_v}, run={run_v}: ok={ok}, msg={msg!r}"


def _in_plane_below_floor_case(tmp_path, comparator=None):
    """Columns 6 and 7 (shear_stress_2 and _3) zero in the reference; the run
    holds 5e-10 in column 6 at row 1. `comparator` is the (gpath, rpath)
    comparison, default the DD station series."""
    def row(t, v):
        r = _station_row(t, -9.0, -9.0)
        r[5], r[6] = v, 0.0
        return r
    ref = [row(1.0, 0.0), row(2.0, 0.0)]
    run = [row(1.0, 0.0), row(2.0, 5e-10)]
    return (comparator or _dd_station_series)(*_station_pair(tmp_path, ref, run))


def _series_floor_case(tmp_path, col, ref_v, run_v, fname="fltst_a.txt"):
    """compare_series (the onfault comparison) on an 11-column pair equal
    except column `col` (0-based) at row 1, built as _floor_case builds it."""
    def row(t, v):
        r = _station_row(t, -9.0, -9.0)
        if col in (2, 6):
            r[col - 1] = 0.0
        r[col] = v
        return r
    g, r = str(tmp_path / "ref" / fname), str(tmp_path / "run" / fname)
    for p, rows in ((g, [row(1.0, ref_v), row(2.0, ref_v)]), (r, [row(1.0, ref_v), row(2.0, run_v)])):
        os.makedirs(os.path.dirname(p))
        _write_station_file(p, rows)
    return C.compare_series(g, r)


@pytest.mark.e2e_fast
@pytest.mark.parametrize("col,name,ref_v,below,above", _FLOOR_CASES)
def test_compare_series_fltst_difference_under_floor_passes(tmp_path, col, name, ref_v, below, above):
    ok, msg = _series_floor_case(tmp_path, col, ref_v, below)
    assert ok, f"{name}: a difference under its floor must pass: {msg}"


@pytest.mark.e2e_fast
@pytest.mark.parametrize("col,name,ref_v,below,above", _FLOOR_CASES)
def test_compare_series_fltst_difference_above_floor_fails(tmp_path, col, name, ref_v, below, above):
    ok, msg = _series_floor_case(tmp_path, col, ref_v, above)
    assert not ok and f"column {col+1}" in msg, f"{name}: ok={ok}, msg={msg!r}"


@pytest.mark.e2e_fast
@pytest.mark.parametrize("ref_v,run_v", [(float("-inf"), -20.0), (float("nan"), -20.0), (-20.0, float("nan"))])
def test_compare_series_fltst_non_finite_log10_slip_rate_is_never_floored(tmp_path, ref_v, run_v):
    ok, msg = _series_floor_case(tmp_path, 4, ref_v, run_v)
    assert not ok and "column 5" in msg, f"ref={ref_v}, run={run_v}: ok={ok}, msg={msg!r}"


@pytest.mark.e2e_fast
@pytest.mark.parametrize("col,name,ref_v,below,above", _FLOOR_CASES)
def test_compare_series_srfst_file_gets_no_floor(tmp_path, col, name, ref_v, below, above):
    ok, msg = _series_floor_case(tmp_path, col, ref_v, below, fname="srfst_a.txt")
    assert not ok and f"column {col+1}" in msg, f"{name}: ok={ok}, msg={msg!r}"


@pytest.mark.e2e_fast
def test_compare_series_fltst_floors_do_not_reach_in_plane_columns(tmp_path):
    """As the DD test of the same name: 0 against 5e-10 in column 6, under the
    shear_stress_3 floor, must fail."""
    ok, msg = _in_plane_below_floor_case(tmp_path, C.compare_series)
    assert not ok and "column 6" in msg, f"ok={ok}, msg={msg!r}"


@pytest.mark.e2e_fast
def test_compare_series_fltst_nine_column_layout_gets_the_floors(tmp_path):
    """9 columns (BP5/BP7): slip_3 is column 3, slip rate along dip (not
    logged) column 5, tstk/tdip columns 6 and 7; shear_stress_3 is column 7."""
    def write(p, v):
        os.makedirs(os.path.dirname(p))
        with open(p, "w") as f:
            for t in (1.0, 2.0):
                row = [t, 1.0, 0.0, 1.0, 0.0, 1.0e6, v if t == 2.0 else 0.0, 1.0, -5.0]
                f.write("".join(f"{x:21.13E}" if k == 0 else f"{x:15.6E}" for k, x in enumerate(row)) + "\n")
    g, r, a = (str(tmp_path / d / "fltst_n.txt") for d in ("ref", "run", "far"))
    write(g, 0.0), write(r, 5e-10), write(a, 3e-9)
    assert C.compare_series(g, r)[0] and not C.compare_series(g, a)[0]


@pytest.mark.e2e_fast
def test_kill_group_and_collect_returns_output_of_a_killed_group():
    p = subprocess.Popen(["bash", "-c", "echo before; sleep 60"], stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT, text=True, start_new_session=True)
    assert p.stdout.readline() == "before\n"
    _kill_group_and_collect(p, 10)
    assert p.returncode == -9, f"returncode {p.returncode}"


@pytest.mark.e2e_fast
def test_kill_group_and_collect_fails_when_a_member_escaped_the_group():
    """The launcher's child is in its own session, so killpg misses it and it
    holds the pipe: a bare communicate() after killpg waits for its 60 s."""
    import signal
    child = ("import subprocess, time; "
             "c = subprocess.Popen(['sleep', '60'], start_new_session=True); "
             "print(c.pid, flush=True); time.sleep(60)")
    p = subprocess.Popen(["python3", "-I", "-c", child], stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT, text=True, start_new_session=True)
    escaped = int(p.stdout.readline())
    try:
        with pytest.raises(pytest.fail.Exception, match=f"process group {p.pid} still holds"):
            _kill_group_and_collect(p, 2)
    finally:
        os.kill(escaped, signal.SIGKILL)
        p.communicate(timeout=10)


@pytest.mark.e2e_fast
@pytest.mark.parametrize("ref_value", [float("inf"), float("-inf"), float("nan")])
def test_dd_profile_series_non_finite_reference_is_a_failure_naming_file_and_column(tmp_path, ref_value):
    """Same hole as the station comparison: an inf in a reference column made
    its scale inf and the column passed whatever the run held (rule 2)."""
    import numpy as np
    g, r = str(tmp_path / "ref" / "p1output.txt"), str(tmp_path / "run" / "p1output.txt")
    for p, tdip in ((g, ref_value), (r, 1.0)):
        os.makedirs(os.path.dirname(p))
        np.savetxt(p, np.array([[1.0, 0.1, 0.0, 5.0, 0.0], [2.0, 0.2, 0.0, 6.0, tdip]]))
    ok, msg = _dd_profile_series(g, r)
    assert not ok and "p1output.txt" in msg and "column 5" in msg, f"ok={ok}, msg={msg!r}"


@pytest.mark.e2e_fast
def test_station_velocity_uses_the_step_taken():
    """Regression (1.24.12): off-fault station velocity is displacement change
    over the step actually taken.

    The solver loops reassigned dtev1 to the NEXT step's size before dividing
    the nodal displacement change by it, so with a variable step every
    free-node velocity was off by (next dt)/(this dt). A constant-dt case such
    as test.bp5.qdc.2000 at its default xi hides it, so xi is raised until the
    step grows. SEAS station files interleave displacement and velocity per
    component after the time column. The history is recorded before the step's
    own displacement is stored, so the displacement of row i is the state at
    the end of step i-1 (time t_{i-1}), and the velocity of row i must equal
    (d_i - d_{i-1}) / (t_{i-1} - t_{i-2}): the change over the step that made
    it. The unfixed solver divides by t_i - t_{i-1}, the next step. Row 1
    divides by the first step, which is t_0 (the step starts at time 0).
    Tolerance per entry: C.station_velocity_ratio (rounding-aware).
    """
    import numpy as np
    run_dir = C.run_case("test.bp5.qdc.2000",
                         {"xi": 0.5, "nstep": 101, "nt_out": 101, "HPC_ncpu": 2},
                         os.path.join(WORK_ROOT, "e2e.station_velocity"))
    files = C.station_files(run_dir, "srfst_strk")
    assert files, "run wrote no off-fault station files"
    info = json.load(open(os.path.join(run_dir, "runInfo.json")))
    print(f"\nbinary version {info['version']} (run_dir {run_dir})")
    worst, worst_norm, n_tight = 0.0, 0.0, 0
    for f in files:
        t = np.loadtxt(os.path.join(run_dir, f))[:, 0]
        assert np.diff(t).max() > 1.5 * np.diff(t).min(), \
            f"{f}: time step never varied; the test cannot see the divisor"
        w, wn, nt = C.station_velocity_ratio(os.path.join(run_dir, f))
        worst, worst_norm, n_tight = max(worst, w), max(worst_norm, wn), n_tight + nt
    print(f"station velocity / finite-difference: max |ratio - 1| = {worst:.3e}, "
          f"max |ratio - 1| / allowed = {worst_norm:.3e}, tight entries = {n_tight}")
    assert n_tight > 0, "no entry is tight enough to discriminate the divisor"
    assert worst_norm <= 1.0, \
        f"station velocity disagrees with d(disp)/dt by {worst:.3e} ({worst_norm:.2f}x allowed)"
