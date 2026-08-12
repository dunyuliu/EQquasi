"""End-to-end: BP5, BP5-dip90 and BP7 station time series against gold.

The existing e2e tier (test_bp5_bp7_regression.py) locks a fault-plane
snapshot at step 101 -- useful, but at 101 steps none of these three
benchmarks has nucleated (see scripts/plotRuptureTime.py's docstring), so a
101-step station series is nearly flat and would lock in almost nothing about
the interseismic-to-nucleation physics stations are meant to exercise.

This runs each test compset's *first full cycle* (Q0: istart=1, to its own
`exit_slip_rate = 1e-3` exit rather than a step cap) and checks scalar
quantities read off the resulting station files against
reference/<bench>/gold/summary.json -- the same "read named numbers off the
raw run, not off a subsampled CSV" pattern test_bp8_against_gold.py uses. The
committed fltst_strk*.csv files are for plotting/eyeballing (there is no
existing multi-station overlay tool for BP5/BP7, unlike BP8's
plotAgainstGold.py compare_bp8 -- see the session report), not for this
numeric check.

The 101-step fault-plane snapshot gold is untouched by this file.

Opt-in, like the rest of e2e:

    python3 -m pytest -m e2e tests/
"""

import glob
import json
import os
import re
import shutil
import subprocess
import time

import numpy as np
import pytest

from conftest import ROOT

BENCH_COMPSET = {
    "bp5": "test.bp5.qdc",
    "bp5.dip90": "test.bp5.qdc.dip90",
    "bp7": "test.bp7.qdc",
}

RUN_TIMEOUT_S = 4800  # 80 min; BP5's Q0 alone measured ~30-45 min on this box

RTOL = 1e-4
ATOL_LOG = 1e-3


def _env():
    env = dict(os.environ)
    env["EQQUASIROOT"] = str(ROOT)
    env["PATH"] = f"{ROOT}/bin:{ROOT}/scripts:" + env.get("PATH", "")
    env["OMP_NUM_THREADS"] = "1"
    return env


def _gold(bench):
    path = os.path.join(str(ROOT), "reference", bench, "gold", "summary.json")
    data = json.load(open(path))
    return {k: v for k, v in data.items() if k.startswith("fltst_strk")}


def _run_q0(bench):
    compset = BENCH_COMPSET[bench]
    exe = os.path.join(str(ROOT), "bin", "eqquasi")
    if not os.path.exists(exe):
        pytest.skip("bin/eqquasi not built; see README 'Rebuilding just the Fortran'")
    gold = _gold(bench)
    if not gold:
        pytest.skip(f"no station gold for {bench} yet")

    env = _env()
    work = os.path.join(str(ROOT), "work", f"e2e.{bench}.q0")
    shutil.rmtree(work, ignore_errors=True)
    subprocess.run(["create.newcase", work, compset], cwd=str(ROOT), env=env,
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    udp = os.path.join(work, "user_defined_params.py")
    s = open(udp).read()
    new_s = re.sub(r"^par\.nstep\s*=.*$",
                   "par.nstep       = 10000 # raised so Q0 can reach its own exit_slip_rate exit",
                   s, count=1, flags=re.MULTILINE)
    assert new_s != s, f"{bench}: could not find par.nstep to override"
    open(udp, "w").write(new_s)

    subprocess.run(["./case.setup"], cwd=work, env=env, check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    run_sh = os.path.join(work, "run.sh")
    txt = re.sub(r"mpirun\.openmpi -np \d+ eqquasi", f"mpirun -np 1 {exe}",
                open(run_sh).read())
    open(run_sh, "w").write(txt)
    open(os.path.join(work, "currentcycle.txt"), "w").write("1\n")

    t0 = time.time()
    with open(os.path.join(work, "r.log"), "w") as out, \
         open(os.path.join(work, "t.log"), "w") as err:
        proc = subprocess.Popen(["bash", "run.sh"], cwd=work, env=env,
                                stdout=out, stderr=err, start_new_session=True)
        while proc.poll() is None:
            if time.time() - t0 > RUN_TIMEOUT_S:
                proc.kill()
                pytest.fail(f"{bench} Q0 exceeded {RUN_TIMEOUT_S} s")
            time.sleep(10)
    assert proc.returncode == 0, f"solver exited {proc.returncode}; see {work}/t.log"
    wall = time.time() - t0
    q0 = os.path.join(work, "Q0")
    assert os.path.isdir(q0), f"no Q0/ produced under {work}"
    return q0, wall


@pytest.fixture(scope="module", params=list(BENCH_COMPSET))
def q0_run(request):
    bench = request.param
    q0, wall = _run_q0(bench)
    return bench, q0, wall


def test_q0_station_matches_gold(q0_run):
    bench, q0, wall = q0_run
    gold = _gold(bench)
    bad = []
    for fname, g in gold.items():
        path = os.path.join(q0, f"{fname}.txt")
        if not os.path.exists(path):
            bad.append(f"{fname}: missing from the run")
            continue
        a = np.loadtxt(path)
        checks = [
            ("t_end_s", a[-1, 0], RTOL, False),
            ("slip_strike_final_m", a[-1, 1], RTOL, False),
            ("slip_rate_strike_peak_ms", a[:, 3].max(), RTOL, False),
            ("state_log10_final", a[-1, 8], ATOL_LOG, True),
        ]
        for key, got, tol, is_abs in checks:
            exp = g[key]
            ok = abs(got - exp) <= tol if is_abs else abs(got - exp) <= tol * max(abs(exp), 1e-30)
            if not ok:
                bad.append(f"{fname}.{key}: got {got:.6g}, gold {exp:.6g}")
    assert not bad, (
        f"{bench} Q0 station series diverged from gold (wall time {wall:.0f} s):\n"
        + "\n".join(bad)
    )
