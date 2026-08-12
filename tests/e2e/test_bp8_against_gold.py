"""End-to-end: build, run BP8-QD-GS at dx = 50 m, and diff against the gold.

This is the tier the suite was missing. Everything else here is static -- it
reads source and compset files and never executes the solver, so it cannot catch
a change that alters the physics while leaving the text intact. Exactly that
happened twice in this project: the pore-pressure boundary treatment and the
initial condition both changed the answer by tens of percent without tripping a
single guard.

The oracle is `reference/bp8/gold/summary.json`, frozen from EQquasi v1.4.7:
Table 1's `tau^0` = 14.6 MPa with `theta_0` derived so the fault starts at
`V_init`, and the finite-volume Omega_f boundary. That result agrees with
`taehoKim_ref` to -2 % on slip at the injection point, +5 % at 200 m, and 0.7 %
on the late-time pore pressure.

Opt-in, because it builds the code and runs for ~20 minutes:

    python3 -m pytest -m e2e tests/            # this tier only
    python3 -m pytest -m "not e2e" tests/      # everything else, the default

Tolerances are tight on purpose. This reruns a deterministic case with the same
binary and mesh, so agreement should be to round-off; anything looser would let
a real physics change through. If a legitimate change moves these numbers,
regenerate the gold deliberately and say so in the commit -- do not widen the
tolerance.
"""

import glob
import json
import os
import shutil
import subprocess
import time

import numpy as np
import pytest

from conftest import ROOT

GOLD = os.path.join(str(ROOT), "reference", "bp8", "gold", "summary.json")
GOLD_DIR = os.path.join(str(ROOT), "reference", "bp8", "gold")
WORK = os.path.join(str(ROOT), "work", "e2e.bp8.gold")

# The nine on-fault stations x2, x3 in {-200, 0, 200} m (section 4.1). All nine
# now have gold: the first four were frozen from v1.4.7's work/bp8.sub147, the
# remaining five were added from a v1.6.0 rerun of the same configuration,
# confirmed bit-identical to the first four before being trusted as gold.
ALL_STATIONS = [
    f"{s:+04d}dp{d:+04d}" for s in (-200, 0, 200) for d in (-200, 0, 200)
]

PROFILE_FILES = [
    f"{q}_{line}.dat"
    for q in ("slip_2", "slip_3", "shear_stress_2", "shear_stress_3", "pore_pressure")
    for line in ("strike", "depth")
]

# Deterministic rerun of the same case: agreement should be near round-off.
RTOL_SLIP = 1e-4
ATOL_LOG = 1e-3     # log10 quantities
RTOL_PRESSURE = 1e-4

RUN_TIMEOUT_S = 3600


def _env():
    env = dict(os.environ)
    env["EQQUASIROOT"] = str(ROOT)
    env["PATH"] = f"{ROOT}/bin:{ROOT}/scripts:" + env.get("PATH", "")
    env["OMP_NUM_THREADS"] = "1"
    return env


def _read(path):
    rows = []
    for line in open(path):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        tok = line.split()
        try:
            float(tok[0])
        except ValueError:
            continue
        rows.append([float(v) for v in tok])
    return np.array(rows)


@pytest.fixture(scope="module")
def run_dir():
    if not os.path.exists(GOLD):
        pytest.skip(f"no gold reference at {GOLD}")
    exe = os.path.join(str(ROOT), "bin", "eqquasi")
    if not os.path.exists(exe):
        pytest.skip("bin/eqquasi not built; see README 'Rebuilding just the Fortran'")

    env = _env()
    shutil.rmtree(WORK, ignore_errors=True)
    subprocess.run(["create.newcase", WORK, "test.bp8.qdc"],
                   cwd=str(ROOT), env=env, check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    # Match the gold: xi = 0.2, run to fluid_tend rather than an nstep cap.
    udp = os.path.join(WORK, "user_defined_params.py")
    s = open(udp).read()
    for old, new in [("par.xi = ", "par.xi = 0.2 #"),
                     ("par.nstep       = ", "par.nstep       = 8000 #"),
                     ("par.nt_out      = ", "par.nt_out      = 8000 #")]:
        s = "\n".join(new + ln.split("=", 1)[1].split("#")[0].strip() * 0
                      if ln.startswith(old) else ln for ln in s.split("\n"))
    open(udp, "w").write(s)

    subprocess.run(["./case.setup"], cwd=WORK, env=env, check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    open(os.path.join(WORK, "currentcycle.txt"), "w").write("1\n")

    with open(os.path.join(WORK, "r.log"), "w") as out, \
         open(os.path.join(WORK, "t.log"), "w") as err:
        proc = subprocess.Popen(["mpirun", "-np", "1", os.path.join(str(ROOT), "bin", "eqquasi")],
                                cwd=WORK, env=env, stdout=out, stderr=err,
                                start_new_session=True)
    t0 = time.time()
    while proc.poll() is None:
        if time.time() - t0 > RUN_TIMEOUT_S:
            proc.kill()
            pytest.fail(f"BP8 e2e run exceeded {RUN_TIMEOUT_S} s")
        time.sleep(10)
    assert proc.returncode == 0, f"solver exited {proc.returncode}; see {WORK}/t.log"
    assert os.path.exists(os.path.join(WORK, "global.dat")), "no global.dat produced"
    return WORK


@pytest.fixture(scope="module")
def gold():
    return json.load(open(GOLD))


@pytest.mark.parametrize("station", ALL_STATIONS)
def test_station_matches_gold(run_dir, gold, station):
    g = gold[station]
    a = _read(os.path.join(run_dir, f"fltst_strk{station}.dat"))

    assert a[0, 5] == pytest.approx(g["tau0_MPa"], abs=1e-4), "tau^0 changed"
    assert 10 ** a[1, 3] == pytest.approx(g["V_step1"], rel=1e-3), \
        "the fault no longer starts at V_init"
    assert a[-1, 1] * 1000 == pytest.approx(g["slip_mm"], rel=RTOL_SLIP), \
        f"slip at {station} moved: {a[-1,1]*1000:.4f} mm vs gold {g['slip_mm']:.4f}"
    assert a[:, 7].max() == pytest.approx(g["p_peak_MPa"], rel=RTOL_PRESSURE), \
        "peak pore pressure moved"
    assert a[-1, 7] == pytest.approx(g["p_late_MPa"], rel=RTOL_PRESSURE), \
        "late-time pore pressure moved -- check the Omega_f cell fractions"
    assert a[0, 10] == pytest.approx(g["state_t0"], abs=ATOL_LOG), \
        "initial state moved -- check theta_0 in the compset"
    assert a[-1, 10] == pytest.approx(g["state_end"], abs=ATOL_LOG), "final state moved"


def test_global_matches_gold(run_dir, gold):
    g = gold["global"]
    a = _read(os.path.join(run_dir, "global.dat"))
    assert a[-1, 0] / 86400.0 == pytest.approx(g["t_end_d"], rel=1e-4), \
        "run no longer reaches fluid_tend"
    assert a[:, 1].max() == pytest.approx(g["peak_Vmax_log10"], abs=ATOL_LOG)
    assert a[a[:, 1].argmax(), 0] / 86400.0 == pytest.approx(g["peak_Vmax_day"], rel=1e-2)
    assert a[:, 2].max() == pytest.approx(g["peak_moment"], rel=1e-3)


# Every gold file is diffed whole, not sampled. The scalar checks above name
# *which* quantity moved and are worth keeping for that, but they read a
# handful of entries out of ~58000 per station: a change confined to the middle
# of a curve, or to a column nothing samples (slip_3, the Darcy velocities),
# passes them untouched. These two compare the full array, which is possible
# only because the gold now stores the run at full time resolution rather than
# decimated to 500 rows.
FULL_DIFF_ATOL = 1e-6


def _gold_csv(name):
    p = os.path.join(GOLD_DIR, name)
    if not os.path.exists(p):
        pytest.skip(f"no gold at {p}")
    return np.genfromtxt(p, delimiter=",", skip_header=1)


@pytest.mark.parametrize("station", ALL_STATIONS)
def test_station_series_matches_gold_in_full(run_dir, station):
    gold = _gold_csv(f"fltst_strk{station}.csv")
    run = _read(os.path.join(run_dir, f"fltst_strk{station}.dat"))
    assert run.shape == gold.shape, \
        f"{station}: shape {run.shape} vs gold {gold.shape}"
    d = np.max(np.abs(run - gold))
    assert d == pytest.approx(0.0, abs=FULL_DIFF_ATOL), (
        f"station {station} diverged from gold over the full series: "
        f"max|diff| = {d:.3e} at column "
        f"{int(np.unravel_index(np.argmax(np.abs(run - gold)), run.shape)[1]) + 1}"
    )


def test_global_series_matches_gold_in_full(run_dir):
    gold = _gold_csv("global.csv")
    run = _read(os.path.join(run_dir, "global.dat"))
    assert run.shape == gold.shape, \
        f"global: shape {run.shape} vs gold {gold.shape}"
    d = np.max(np.abs(run - gold))
    assert d == pytest.approx(0.0, abs=FULL_DIFF_ATOL), \
        f"global.dat diverged from gold over the full series: max|diff| = {d:.3e}"


def test_fault_snapshot_matches_gold(run_dir):
    """BP5/BP7 have long compared a frozen fault.*.nc; BP8 never had one.

    Confirmed before freezing: work/bp8.sub147 (v1.4.7) and a fresh v1.6.0 rerun
    at this same configuration agree on every variable to max|diff| = 0.0, the
    same bar BP5/BP7 hold their snapshots to (reference/bp5/gold/summary.json).
    A single MPI rank, deterministic solver and fixed seed make that the right
    tolerance -- anything looser would let a real change through silently.
    """
    netCDF4 = pytest.importorskip("netCDF4")

    hits = glob.glob(os.path.join(run_dir, "fault.?????.nc"))
    hits = [h for h in hits if not h.endswith("fault.r.nc")
            and not os.path.basename(h) == "fault.00001.nc"]
    assert hits, f"no fault.*.nc snapshot under {run_dir}"
    run_path = max(hits, key=os.path.getmtime)

    gold_path = os.path.join(GOLD_DIR, "fault.05301.nc")
    if not os.path.exists(gold_path):
        pytest.skip("no fault.05301.nc gold")

    g = netCDF4.Dataset(gold_path)
    r = netCDF4.Dataset(run_path)
    bad = []
    for v in g.variables:
        if v.startswith("nid_"):
            continue
        gv = np.asarray(g.variables[v][:]).squeeze()
        rv = np.asarray(r.variables[v][:]).squeeze()
        if gv.shape != rv.shape:
            bad.append(f"{v}: shape {rv.shape} vs gold {gv.shape}")
            continue
        d = np.max(np.abs(gv - rv))
        if d > 0.0:
            bad.append(f"{v}: max|diff| = {d:.3e}")
    assert not bad, "fault-plane snapshot diverged from gold:\n" + "\n".join(bad)


@pytest.mark.parametrize("fname", PROFILE_FILES)
def test_section_43_profile_matches_gold(run_dir, fname):
    """Section 4.3 profiles at the gold's native dx = 50 m (17 nodes), not the
    submission's resampled 81-node grid -- that resampling is
    resampleBP8Profiles.py's job and is covered by test_submission_validator_passes.
    """
    gold_path = os.path.join(GOLD_DIR, fname.replace(".dat", ".csv"))
    if not os.path.exists(gold_path):
        pytest.skip(f"no gold for {fname}")
    gold = np.genfromtxt(gold_path, delimiter=",", skip_header=1)

    run_path = os.path.join(run_dir, fname)
    assert os.path.exists(run_path), f"{fname} not produced by the run"
    run = _read(run_path)

    assert run.shape == gold.shape, (
        f"{fname}: shape {run.shape} vs gold {gold.shape}"
    )
    # Row 0 is [0, 0, coords...]; later rows are [t, max_slip_rate, values...].
    # Deterministic rerun of the same config: hold to the same bar as the
    # station and global checks above.
    assert np.max(np.abs(run - gold)) == pytest.approx(0.0, abs=1e-6), (
        f"{fname} diverged from gold: max|diff| = {np.max(np.abs(run - gold)):.3e}"
    )


def test_submission_validator_passes(run_dir):
    """The run must also be uploadable, not merely numerically right."""
    out = os.path.join(str(ROOT), "work", "e2e.sub")
    shutil.rmtree(out, ignore_errors=True)
    subprocess.run(["python3", os.path.join(str(ROOT), "scripts", "resampleBP8Profiles.py"),
                    run_dir, out], check=True, stdout=subprocess.DEVNULL)
    r = subprocess.run(["python3", os.path.join(str(ROOT), "scripts", "checkBP8Submission"), out],
                       capture_output=True, text=True)
    assert r.returncode == 0, f"checkBP8Submission rejected the e2e output:\n{r.stdout}"
