"""End-to-end: the two-fault step-over (case_input/test.stepover.qdc), the
suite's only ntotft > 1 gate.

Every other e2e gate -- BP5, BP5-dip90, BP7, BP8 -- runs a single-fault
compset. Three multi-fault bugs reached master and every one was found only
by a person running a new case by hand, never by the suite:

  1. accumulator aliasing reaching dtev (state advanced with the wrong step
     on one fault)
  2. netcdf_write_on_fault reading 735 slots of uninitialised memory per file
     for a fault the writer did not expect
  3. a fault plane falling between mesh lines: fault A got zero nodes, and
     the run still completed 4240 steps and reported 11 m of slip on a
     phantom single fault, writing an all-zero slab for the missing one

None of these looked like a crash. All three produced a run that completed
and looked plausible from the outside. This gate is built specifically to
make each one visible: a missing/empty/constant/NaN slab is asserted on
directly, not inferred from a scalar summary.

The compset itself encodes the check for bug #1's sibling class (a swapped or
mis-indexed per-fault read): fault 0's init_norm is -25 MPa, fault 1's is
-27 MPa, chosen to be distinguishable but both physically unremarkable.

Opt-in, like the rest of e2e:

    python3 -m pytest -m e2e tests/
"""

import glob
import json
import os
import shutil
import subprocess

import numpy as np
import pytest

from conftest import ROOT

COMPSET = "test.stepover.qdc"
GOLD_DIR = os.path.join(str(ROOT), "reference", "stepover")
WORK = os.path.join(str(ROOT), "work", "e2e.stepover")
SNAPSHOT = "fault.00101.nc"

RUN_TIMEOUT_S = 1800


def _env():
    env = dict(os.environ)
    env["EQQUASIROOT"] = str(ROOT)
    env["PATH"] = f"{ROOT}/bin:{ROOT}/scripts:" + env.get("PATH", "")
    env["OMP_NUM_THREADS"] = "1"
    return env


@pytest.fixture(scope="module")
def run_dir():
    gold_path = os.path.join(GOLD_DIR, SNAPSHOT)
    if not os.path.exists(gold_path):
        pytest.skip(f"no gold at {gold_path}")
    exe = os.path.join(str(ROOT), "bin", "eqquasi")
    if not os.path.exists(exe):
        pytest.skip("bin/eqquasi not built")

    env = _env()
    shutil.rmtree(WORK, ignore_errors=True)
    subprocess.run(["create.newcase", WORK, COMPSET], cwd=str(ROOT), env=env,
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    subprocess.run(["./case.setup"], cwd=WORK, env=env, check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    run_sh = os.path.join(WORK, "run.sh")
    import re
    txt = re.sub(r"mpirun\.openmpi -np \d+ eqquasi", f"mpirun -np 1 {exe}",
                open(run_sh).read())
    open(run_sh, "w").write(txt)
    open(os.path.join(WORK, "currentcycle.txt"), "w").write("1\n")

    r = subprocess.run(["bash", "run.sh"], cwd=WORK, env=env,
                       capture_output=True, text=True, timeout=RUN_TIMEOUT_S)
    assert r.returncode == 0, f"solver exited {r.returncode}:\n{r.stdout[-3000:]}\n{r.stderr[-2000:]}"
    # run.sh's cycle loop moves every output file into ./Q<cycle-1>/;
    # istart = 1 means cycle 1, i.e. Q0 (same convention as the BP5/BP7/BP8
    # e2e fixtures' own run.sh-driven runs).
    q0 = os.path.join(WORK, "Q0")
    hits = glob.glob(os.path.join(q0, SNAPSHOT))
    assert hits, f"no {SNAPSHOT} produced under {q0}"
    return q0


def _load(path):
    netCDF4 = pytest.importorskip("netCDF4")
    d = netCDF4.Dataset(path)
    return {v: np.asarray(d.variables[v][:]) for v in d.variables}


VARIABLES = ["shear_strike", "shear_dip", "effective_normal", "slip_rate",
            "state_variable", "state_normal", "vxm", "vym", "vzm",
            "vxs", "vys", "vzs", "slips", "slipd", "slipn"]


def test_run_writes_a_slab_for_every_declared_fault(run_dir):
    """Bug #3's signature: ntotft = 2 but the run silently wrote for one.

    netcdf_write_on_fault sizes its leading dimension from ntotft, so
    "wrote fewer slabs than declared" is the direct symptom of a fault that
    meshed zero nodes and never should have been allowed to run.
    """
    run = _load(os.path.join(run_dir, SNAPSHOT))
    assert run["nid_fault"].shape[0] == 2, (
        f"expected 2 fault slabs (ntotft=2), got {run['nid_fault'].shape[0]}"
    )
    for v in VARIABLES:
        assert run[v].shape[0] == 2, f"{v}: {run[v].shape[0]} slabs, expected 2"


@pytest.mark.parametrize("ift", [0, 1])
def test_no_slab_is_all_zero_or_contains_nan(run_dir, ift):
    """Bug #2's signature: uninitialised memory reads as either garbage NaN
    or, on a stack that happened to be pre-zeroed, a silently all-zero slab
    that looks like 'nothing happened here yet' rather than 'this never
    initialised'.
    """
    run = _load(os.path.join(run_dir, SNAPSHOT))
    bad = []
    for v in VARIABLES:
        arr = run[v][ift]
        if np.isnan(arr).any():
            bad.append(f"{v}[fault {ift}]: contains NaN")
        elif not np.any(arr):
            bad.append(f"{v}[fault {ift}]: all zero")
    assert not bad, "\n".join(bad)


def test_the_two_faults_have_distinct_effective_normal_stress(run_dir):
    """The compset sets init_norm to -25 MPa (fault 0) and -27 MPa (fault 1)
    specifically so a swapped or aliased per-fault read is visible as a
    value on the wrong fault, not merely a value that looks plausible either
    way (reference/stepover/summary.json's per_fault_design).
    """
    run = _load(os.path.join(run_dir, SNAPSHOT))
    en = run["effective_normal"]
    m0, m1 = en[0].mean(), en[1].mean()
    assert m0 == pytest.approx(-25.0e6, rel=0.05), f"fault 0 mean effective_normal = {m0:.3e}"
    assert m1 == pytest.approx(-27.0e6, rel=0.05), f"fault 1 mean effective_normal = {m1:.3e}"
    assert abs(m0 - m1) > 1.0e6, (
        "fault 0 and fault 1 effective_normal are not distinguishable -- "
        "a swap between the two faults would not be caught"
    )


def test_field_snapshot_matches_gold_per_fault(run_dir):
    """The regression lock proper: every variable, every fault, max|diff| = 0.

    Mirrors BP5/BP7's tolerance (reference/bp5/summary.json). A fault
    present in gold but missing from the run (bug #3's exact shape) is named,
    not merged away into an aggregate that could pass by averaging it out.
    """
    gold = _load(os.path.join(GOLD_DIR, SNAPSHOT))
    run = _load(os.path.join(run_dir, SNAPSHOT))
    bad = []
    for v in VARIABLES:
        if v not in run:
            bad.append(f"{v}: missing from the run entirely")
            continue
        g, r = gold[v], run[v]
        if g.shape != r.shape:
            bad.append(f"{v}: shape {r.shape} vs gold {g.shape}")
            continue
        for ift in range(g.shape[0]):
            d = np.max(np.abs(g[ift] - r[ift]))
            if d > 0.0:
                bad.append(f"{v}[fault {ift}]: max|diff| = {d:.3e}")
    assert not bad, "step-over snapshot diverged from gold:\n" + "\n".join(bad)
