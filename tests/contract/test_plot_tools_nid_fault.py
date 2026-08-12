"""Guard: the plotting tools survive the 3-D on-fault layout, per fault.

`netcdf_write_on_fault` (src/netcdf_io.f90) stores every on-fault variable as
(nid_fault, nid_dip, nid_strike). That shape has now broken consumers three
separate times (check.test.py gold, plotAgainstGold.py's key-collision, the
plotting utilities), each time because a consumer assumed 2-D or silently
reduced the fault dimension -- xarray's .plot() on a 3-D array happily draws a
histogram instead of a map. This test runs every converted tool, as a user
would, against a synthetic two-fault case and asserts:

  * exit code 0 and a figure per fault (multi-fault output is never dropped),
  * the legacy 2-D layout still loads (pre-mfault files remain readable),
  * a directory with none of the required files produces a named error, not a
    traceback,
  * `-h` documents the process-every-cycle default (that behaviour is not
    discoverable any other way).

The tools are exercised through subprocess from inside the case directory with
no arguments -- the primary workflow -- so cycle discovery (Q*), parameter
resolution (user_defined_params.py from the cycle's parent) and the nid_fault
handling are all on the hook together.
"""

import os
import subprocess
import sys
import textwrap

import numpy as np
import pytest

from tests.conftest import ROOT

SCRIPTS = ROOT / "scripts"
ND, NS, NF = 5, 7, 2          # dip nodes, strike nodes, faults
NSTEPS = 101


def write_fault_nc(path, legacy_2d=False, slip=0.0):
    nc = pytest.importorskip("netCDF4")
    d = nc.Dataset(path, "w")
    d.createDimension("nid_dip", ND)
    d.createDimension("nid_strike", NS)
    dims = ("nid_dip", "nid_strike")
    if not legacy_2d:
        d.createDimension("nid_fault", NF)
        dims = ("nid_fault",) + dims
    for name, val in [("shear_strike", 15e6), ("shear_dip", 0.0),
                      ("effective_normal", -25e6), ("slip_rate", 1e-9),
                      ("state_variable", 1e8), ("state_normal", 25e6),
                      ("slips", slip), ("slipd", 0.0)]:
        v = d.createVariable(name, "f8", dims)
        v[:] = val
    d.close()


def make_case(root, legacy_2d=False):
    """A minimal two-fault case laid out the way run.sh leaves it: parameters
    in the case root, results in Q0/."""
    case = root / ("case2d" if legacy_2d else "case")
    q0 = case / "Q0"
    q0.mkdir(parents=True)
    (case / "user_defined_params.py").write_text(textwrap.dedent(f"""\
        import numpy as np
        class _P: pass
        par = _P()
        par.casename = "synthetic"
        par.bp = 5
        par.ntotft = {NF}
        par.init_norm = -25.0e6
        par.fric_rsf_r0 = 0.6
        par.fxmin, par.fxmax = -6.0e3, 6.0e3
        par.fzmin, par.fzmax = -8.0e3, 0.0
        par.dx = par.dz = 2.0e3
        par.fx = np.linspace(par.fxmin, par.fxmax, {NS})
        par.fz = np.linspace(par.fzmin, par.fzmax, {ND})
        par.on_fault_vars = np.zeros(({ND}, {NS}, 100))
        """))
    for step, slip in ((1, 0.0), (NSTEPS, 1.0)):
        write_fault_nc(q0 / f"fault.{step:05d}.nc", legacy_2d, slip)
    write_fault_nc(q0 / "fault.r.nc", legacy_2d, 1.0)
    t = np.linspace(0.065, 300.0, NSTEPS)
    g = np.zeros((NSTEPS, 7))
    g[:, 0], g[:, 1] = t, 1e-9
    np.savetxt(q0 / "global.dat", g)
    (q0 / "tdyna.txt").write_text("0.065\n")
    # cplot: NF faults' nodes stacked, 16 columns, fnft partially ruptured.
    xs = np.linspace(-6e3, 6e3, NS)
    zs = np.linspace(-8e3, 0.0, ND)
    rows = []
    for ift in range(NF):
        for z in zs:
            for x in xs:
                fnft = 2.5 if (ift == 0 and x < 0) else -1000.0
                rows.append([x + ift * 20e3, z] + [0.0] * 13 + [fnft])
    np.savetxt(q0 / "cplot_EQquasi.txt", np.array(rows))
    return case


def run_tool(name, case, *args):
    r = subprocess.run(
        [sys.executable, str(SCRIPTS / name), *args],
        cwd=str(case), capture_output=True, text=True, timeout=300,
        env={**os.environ, "MPLBACKEND": "Agg"})
    assert "Traceback" not in r.stderr, \
        f"{name} crashed with a traceback:\n{r.stderr}"
    return r


@pytest.fixture(scope="module")
def case(tmp_path_factory):
    return make_case(tmp_path_factory.mktemp("plotcase"))


def test_plotOnFaultVars_plots_every_fault(case):
    r = run_tool("plotOnFaultVars", case)
    assert r.returncode == 0, r.stderr or r.stdout
    for ift in range(NF):
        assert (case / "Q0" / f"fault.{NSTEPS:05d}.nc.f{ift}.png").exists(), \
            f"fault {ift} was dropped from a {NF}-fault snapshot"


def test_plotProfiles_plots_every_fault(case):
    r = run_tool("plotProfiles", case)
    assert r.returncode == 0, r.stderr or r.stdout
    for ift in range(NF):
        assert (case / "Q0" / f"cVerticalProfile.slips.f{ift}.png").exists()


def test_plotProfile_runs(case):
    r = run_tool("plotProfile.py", case)
    assert r.returncode == 0, r.stderr or r.stdout
    assert (case / "Q0" / "profile.slips.horizontal.png").exists()


def test_plotAccumulated_runs(case):
    r = run_tool("plotAccumulated", case)
    assert r.returncode == 0, r.stderr or r.stdout
    assert (case / "Q0" / "accumulatedSlip.horizontal.png").exists()


def test_plotRuptureTime_runs(case):
    r = run_tool("plotRuptureTime.py", case)
    assert r.returncode == 0, r.stderr or r.stdout
    assert (case / "Q0" / "rupture_time.png").exists()
    assert "ruptured" in r.stdout


def test_plotPeakSliprateTime_runs(case):
    r = run_tool("plotPeakSliprateTime.py", case)
    assert r.returncode == 0, r.stderr or r.stdout
    assert (case / "peak_slip_rate_vs_time.png").exists()


def test_out_of_range_fault_index_is_a_named_error(case):
    r = run_tool("plotProfile.py", case, "--fault", str(NF))
    assert r.returncode != 0
    assert "out of range" in (r.stderr + r.stdout)


def test_legacy_2d_snapshots_still_load(tmp_path):
    case = make_case(tmp_path, legacy_2d=True)
    r = run_tool("plotOnFaultVars", case)
    assert r.returncode == 0, r.stderr or r.stdout
    assert (case / "Q0" / f"fault.{NSTEPS:05d}.nc.png").exists()


def test_empty_directory_is_a_named_error(tmp_path):
    (tmp_path / "nothing").mkdir()
    r = run_tool("plotRuptureTime.py", tmp_path / "nothing")
    assert r.returncode != 0
    assert "required files" in (r.stderr + r.stdout)


def test_help_documents_cycle_discovery(case):
    for tool in ("plotOnFaultVars", "plotRuptureTime.py", "plotProfile.py",
                 "plotProfiles", "plotAccumulated", "plotPeakSliprateTime.py",
                 "plotCompareStressRoughness"):
        r = run_tool(tool, case, "-h")
        assert r.returncode == 0
        assert "all cycles" in r.stdout, \
            f"{tool} -h does not explain the process-every-cycle default"
