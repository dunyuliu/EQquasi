"""Every benchmark's gold must carry both oracle shapes: a fault-plane
snapshot and a station time series (plus, for BP8, the ten section-4.3
profiles). Before this session BP5/BP7/BP5-dip90 had only the snapshot and
BP8 had only stations; a reference file nothing reads is dead weight, so
this pins what must exist rather than trusting the directory listing.

BP8 got both shapes plus the profiles this session (see reference/test.bp8.qdc.gs.10/
summary.json's "extended_2026-08-12" block for provenance). BP5's Q0 -- the
first full seismic cycle, needed for a station series that actually spans
nucleation rather than 101 steps of flat interseismic loading -- was run and
observed through peak slip rate (~1.14 m/s) and well into decay, but did not
reach its exit_slip_rate = 1e-3 exit inside this session under a heavily
loaded box (load average 20-40 on 64 cores); BP5-dip90 and BP7 were not
started at all, since only one such run fits on the box at a time. Station
gold for all three is deferred -- see the session report -- rather than
frozen from a truncated run, which would silently lock in an arbitrary
mid-decay slip-rate snapshot as if it were the cycle's answer.
"""

import glob
import json
import os

import pytest

from conftest import ROOT

def _reference_runs():
    """Discovered, not listed -- see testsys/unit/test_physical_invariants.py."""
    import os
    out = []
    base = ROOT / "reference"
    for bench in sorted(p for p in base.iterdir() if p.is_dir()):
        for cand in [bench] + sorted(d for d in bench.iterdir()
                                     if d.is_dir() and d.name not in ("plots", "archive")):
            if any(f.name.startswith("fault.") and f.name.endswith(".nc")
                   for f in cand.iterdir() if f.is_file()):
                out.append(str(cand.relative_to(base)))
    return out


FIELD_BENCHMARKS = tuple(_reference_runs())

BP8_STATIONS = [f"{s:+04d}dp{d:+04d}" for s in (-200, 0, 200) for d in (-200, 0, 200)]

BP8_PROFILES = [
    f"{q}_{line}"
    for q in ("slip_2", "slip_3", "shear_stress_2", "shear_stress_3", "pore_pressure")
    for line in ("strike", "depth")
]


def gold_dir(bench):
    return ROOT / "reference" / bench


@pytest.mark.parametrize("bench", FIELD_BENCHMARKS)
def test_field_benchmark_still_has_its_snapshot(bench):
    """Every reference must carry at least one fault-plane snapshot.

    Not `fault.00101.nc` specifically: a reference is whatever its run
    produced, and the step number depends on nt_out and how long the cycle
    ran. bp5/cycle0 ends at 04483, bp8 at 05301.
    """
    d = gold_dir(bench)
    snaps = sorted(d.glob("fault.*.nc"))
    assert snaps, f"{bench} has no fault-plane snapshot"
    assert any(not f.name.endswith(".r.nc") for f in snaps), \
        f"{bench} has only a restart file, no numbered snapshot"


@pytest.mark.parametrize("station", BP8_STATIONS)
def test_bp8_has_all_nine_stations(station):
    d = gold_dir("test.bp8.qdc.gs.10")
    assert (d / f"fltst_strk{station}.csv").is_file()
    summary = json.loads((d / "summary.json").read_text())
    assert station in summary, f"summary.json carries no numbers for station {station}"


@pytest.mark.parametrize("profile", BP8_PROFILES)
def test_bp8_has_all_ten_section43_profiles(profile):
    d = gold_dir("test.bp8.qdc.gs.10")
    assert (d / f"{profile}.csv").is_file(), (
        f"reference/test.bp8.qdc.gs.10/{profile}.csv is missing; "
        "resampleBP8Profiles.py and checkBP8Submission both exercise this "
        "file with no oracle to catch a regression in it"
    )


def test_bp8_has_a_fault_plane_snapshot():
    """BP5/BP7 have compared a fault-plane snapshot for a long time; BP8 never
    did. Filename differs from BP5/BP7 (step 05301, not 00101) because BP8 is
    a single aseismic run, not a multi-cycle one -- see reference/test.bp8.qdc.gs.10/README.md.
    """
    d = gold_dir("test.bp8.qdc.gs.10")
    hits = glob.glob(str(d / "fault.*.nc"))
    assert hits, "reference/test.bp8.qdc.gs.10 has no fault-plane netCDF snapshot"
    csv_hits = glob.glob(str(d / "fault.*.csv"))
    assert csv_hits, "reference/test.bp8.qdc.gs.10 has no flattened CSV of the fault-plane snapshot"


def test_bp8_fault_snapshot_csv_matches_the_netcdf():
    """reference/test.bp8.qdc.gs.10/fault.05301.csv is the flattened, human-eyeballable
    twin of fault.05301.nc (same convention as BP5/BP7's fault.00101.csv). A
    reader that only checks the .nc would miss the .csv going stale if either
    is regenerated without the other -- this cross-checks them row for row.
    """
    netCDF4 = pytest.importorskip("netCDF4")
    import numpy as np

    d = gold_dir("test.bp8.qdc.gs.10")
    nc_path = d / "fault.05301.nc"
    csv_path = d / "fault.05301.csv"
    assert csv_path.is_file(), "reference/test.bp8.qdc.gs.10/fault.05301.csv is missing"

    ds = netCDF4.Dataset(nc_path)
    names = [v for v in ds.variables if not v.startswith("nid_fault")]
    dip_n = ds.variables["nid_dip"].shape[0]
    strike_n = ds.variables["nid_strike"].shape[0]

    header = open(csv_path).readline().strip().split(",")
    assert header == names, (
        f"fault.05301.csv column order {header} does not match the netCDF "
        f"variable order {names}"
    )

    rows = np.genfromtxt(csv_path, delimiter=",", skip_header=1)
    assert rows.shape[0] == dip_n * strike_n, (
        f"fault.05301.csv has {rows.shape[0]} rows, expected "
        f"{dip_n} x {strike_n} = {dip_n * strike_n}"
    )
    for col, v in enumerate(names):
        arr = np.asarray(ds.variables[v][:]).squeeze()
        flat = arr.reshape(-1) if arr.ndim == 2 else np.repeat(arr, strike_n) \
            if v == "nid_dip" else np.tile(arr, dip_n)
        assert np.allclose(rows[:, col], flat), (
            f"fault.05301.csv column '{v}' disagrees with the netCDF"
        )


def test_bp8_global_csv_matches_summary():
    """global.csv's last row must agree with the scalar numbers
    summary.json's "global" block and the e2e test both hold it to."""
    d = gold_dir("test.bp8.qdc.gs.10")
    csv_path = d / "global.csv"
    assert csv_path.is_file(), "reference/test.bp8.qdc.gs.10/global.csv is missing"

    import numpy as np
    rows = np.genfromtxt(csv_path, delimiter=",", skip_header=1)
    g = json.loads((d / "summary.json").read_text())["global"]
    assert rows[-1, 0] / 86400.0 == pytest.approx(g["t_end_d"], rel=1e-3)
    assert rows[:, 1].max() == pytest.approx(g["peak_Vmax_log10"], abs=1e-3)


# Row 6 (2026-09-25): peak_sliprate_per_fault.dat is time, then one peak-V
# column per fault. It is additive to global.dat (whose column 2 stays the
# all-fault max), so a reference that carries it must agree with its own
# global.dat: same rows, same time column, and the max over the fault columns
# is column 2 -- both are maxval() of the same sliprate_arr in faulting.f90 at
# the same step, so within one run they are bit-identical. Across the files of
# a reference they are not: rule 8 keeps the pre-existing global.dat (knox,
# 1.16.0 for the stepover), and the per-fault file was blessed additively from
# a later run, so the two differ by run-to-run MPI reduction noise (~1e-15
# relative, measured 2026-09-25). The gate's own series tolerance (1e-9,
# testsys/e2e/cases.py) is used here for the same reason it is used there; a
# routing bug puts a wrong FAULT's value in a column, which is not a 1e-15
# effect. This is also the rule-8a reader for the file in the fast tier;
# testsys/e2e/cases.py compares it against a fresh run in the e2e tier.
PERFAULT_REFS = [b for b in FIELD_BENCHMARKS
                 if (gold_dir(b) / "peak_sliprate_per_fault.dat").is_file()]


def test_some_reference_carries_peak_sliprate_per_fault():
    assert PERFAULT_REFS, ("no reference holds peak_sliprate_per_fault.dat; the "
                           "per-fault test below would be vacuous")


@pytest.mark.parametrize("bench", PERFAULT_REFS)
def test_peak_sliprate_per_fault_matches_global(bench):
    import sys
    import numpy as np
    sys.path.insert(0, str(ROOT / "script"))
    from seasio import read_array
    d = gold_dir(bench)
    p = np.atleast_2d(read_array(d / "peak_sliprate_per_fault.dat"))
    g = np.atleast_2d(read_array(d / "global.dat"))
    assert p.shape[1] >= 2, f"{bench}: fewer than one fault column"
    assert p.shape[0] == g.shape[0], \
        f"{bench}: {p.shape[0]} rows vs {g.shape[0]} in global.dat"
    assert np.array_equal(p[:, 0], g[:, 0]), f"{bench}: time columns differ"
    vmax = p[:, 1:].max(axis=1)
    assert np.allclose(vmax, g[:, 1], rtol=1e-9, atol=1e-18), \
        (f"{bench}: max over faults != global.dat column 2; worst relative "
         f"difference {np.max(np.abs(vmax - g[:, 1]) / np.maximum(np.abs(g[:, 1]), 1e-30)):.2e}")
    # The compset says how many faults there are; the file must have one
    # column per fault, or a fault is silently missing from the figure.
    import re
    udp = (d / "user_defined_params.py").read_text() if (d / "user_defined_params.py").is_file() else ""
    m = re.search(r"par\.ntotft\s*=\s*(\d+)", udp)
    if m:
        assert p.shape[1] - 1 == int(m.group(1)), \
            f"{bench}: {p.shape[1] - 1} fault columns, par.ntotft = {m.group(1)}"
