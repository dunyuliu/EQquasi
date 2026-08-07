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

import json
import os
import shutil
import subprocess
import time

import numpy as np
import pytest

from conftest import ROOT

GOLD = os.path.join(str(ROOT), "reference", "bp8", "gold", "summary.json")
WORK = os.path.join(str(ROOT), "work", "e2e.bp8.gold")

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


@pytest.mark.parametrize("station", ["+000dp+000", "-200dp+000"])
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


def test_submission_validator_passes(run_dir):
    """The run must also be uploadable, not merely numerically right."""
    out = os.path.join(str(ROOT), "work", "e2e.sub")
    shutil.rmtree(out, ignore_errors=True)
    subprocess.run(["python3", os.path.join(str(ROOT), "scripts", "resampleBP8Profiles.py"),
                    run_dir, out], check=True, stdout=subprocess.DEVNULL)
    r = subprocess.run(["python3", os.path.join(str(ROOT), "scripts", "checkBP8Submission"), out],
                       capture_output=True, text=True)
    assert r.returncode == 0, f"checkBP8Submission rejected the e2e output:\n{r.stdout}"
