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

import json
import os
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
    ("bp5",       "test.bp5.qdc",       {}, "cycle0-step101-fast", "fast"),
    ("bp8",       "test.bp8.qdc",
     {"xi": 0.2, "nstep": 8000, "nt_out": 8000}, "",               "fast"),

    # full: -m e2e only.
    ("bp5",       "test.bp5.qdc",
     {"nstep": 100000, "nt_out": 1000},          "cycle0",         "full"),

    # Recreated later, through this mechanism, from a clean run: BP5-dip90,
    # BP7 and the step-over. Their old references were produced ad hoc, before
    # the workflow rule and before the gold/ flatten, so they are not carried
    # forward -- a reference nobody can regenerate is worse than none. Add a
    # row here when its reference exists and has been verified.
]

RUN_TIMEOUT_S = 7200


def case_id(c):
    name, _, over, ref, tier = c
    return f"{name}.{ref}" if ref else name


def reference_dir(name, ref):
    d = os.path.join(str(ROOT), "reference", name)
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


def run_case(compset, over, workdir):
    """Drive one case through the workflow and return its cycle0 directory.

    create.newcase -> case.setup -> run.sh, never the binary directly: run.sh
    owns the cycle loop, and bypassing it produces output shaped unlike any
    other run's (rule 19).
    """
    e = env()
    exe = os.path.join(str(ROOT), "bin", "eqquasi")
    if not os.path.exists(exe):
        # Not a skip. A skipped benchmark takes the entire e2e gate with it and
        # still reports green -- CI did exactly that: 2 passed, 18 skipped in
        # 8 seconds, having run no benchmark at all. Build first.
        pytest.fail(f"{exe} does not exist. Build before running the e2e "
                    "tier: EQQUASIROOT=$(pwd) MACHINE=<host> make -C src && "
                    "mv src/eqquasi bin/  (or run install.eqquasi.sh)")

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

    shutil.rmtree(workdir, ignore_errors=True)
    step(["create.newcase", workdir, compset], str(ROOT))
    apply_overrides(os.path.join(workdir, "user_defined_params.py"), over)
    step(["./case.setup"], workdir)

    with open(os.path.join(workdir, "run.log"), "w") as log:
        proc = subprocess.Popen(["bash", "run.sh"], cwd=workdir, env=e,
                                stdout=log, stderr=subprocess.STDOUT,
                                start_new_session=True)
    t0 = time.time()
    while proc.poll() is None:
        if time.time() - t0 > RUN_TIMEOUT_S:
            proc.kill()
            pytest.fail(f"{compset} exceeded {RUN_TIMEOUT_S} s")
        time.sleep(10)
    assert proc.returncode == 0, (
        f"run.sh exited {proc.returncode} for {compset}; see {workdir}/run.log")

    cyc = os.path.join(workdir, "cycle0")
    assert os.path.isdir(cyc), f"run.sh produced no cycle0/ for {compset}"
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
    for ext in (".dat", ".txt", ".csv"):
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
    "onfault":  (("fltst_strk*",),                   "series"),
    "offfault": (("srfst_strk*",),                   "series"),
    "global":   (("global.dat", "global.csv"),       "series"),
    "profile":  (("*_strike.dat", "*_strike.csv",
                  "*_depth.dat",  "*_depth.csv"),    "series"),
    "cplot":    (("cplot_EQquasi.*",),               "series"),
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


def compare_series(gpath, rpath, atol=1e-6):
    """(ok, message) for two numeric tables. Full array, not sampled."""
    import numpy as _np
    g, r = read_any(gpath), read_any(rpath)
    if g.shape != r.shape:
        return False, f"shape {r.shape} vs reference {g.shape}"
    d = float(_np.max(_np.abs(g - r)))
    if d > atol:
        col = int(_np.unravel_index(_np.argmax(_np.abs(g - r)), g.shape)[1]) + 1
        return False, f"max|diff| = {d:.3e} at column {col}"
    return True, f"{g.shape[0]} rows x {g.shape[1]} cols, max|diff| = {d:.1e}"


def compare_netcdf(gpath, rpath):
    """(ok, message) per fault, so a missing or aliased fault is named."""
    import numpy as _np
    import netCDF4 as _nc
    g, r = _nc.Dataset(gpath), _nc.Dataset(rpath)
    gd = {k: len(v) for k, v in g.dimensions.items()}
    rd = {k: len(v) for k, v in r.dimensions.items()}
    if gd != rd:
        return False, f"dimensions {rd} vs reference {gd}"
    worst, where = 0.0, ""
    for v in g.variables:
        if v.startswith("nid_"):
            continue
        a, b = _np.asarray(g[v][:], float), _np.asarray(r[v][:], float)
        nf = a.shape[0] if a.ndim == 3 else 1
        for ift in range(nf):
            x = a[ift] if a.ndim == 3 else a
            y = b[ift] if b.ndim == 3 else b
            d = float(_np.max(_np.abs(x - y)))
            if d > worst:
                worst, where = d, f"{v}{f'[fault {ift}]' if nf > 1 else ''}"
    if worst > 0.0:
        return False, f"max|diff| = {worst:.3e} in {where}"
    return True, f"{len(gd)} dims, {len(g.variables)} variables, max|diff| = 0"


COMPARATORS = {"series": compare_series, "netcdf": compare_netcdf}
