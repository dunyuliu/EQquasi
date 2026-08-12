"""Every benchmark's gold must carry both oracle shapes: a fault-plane
snapshot and a station time series (plus, for BP8, the ten section-4.3
profiles). Before this session BP5/BP7/BP5-dip90 had only the snapshot and
BP8 had only stations; a reference file nothing reads is dead weight, so
this pins what must exist rather than trusting the directory listing.

BP8 got both shapes plus the profiles this session (see reference/bp8/gold/
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

FIELD_BENCHMARKS = ("bp5", "bp5.dip90", "bp7")

BP8_STATIONS = [f"{s:+04d}dp{d:+04d}" for s in (-200, 0, 200) for d in (-200, 0, 200)]

BP8_PROFILES = [
    f"{q}_{line}"
    for q in ("slip_2", "slip_3", "shear_stress_2", "shear_stress_3", "pore_pressure")
    for line in ("strike", "depth")
]


def gold_dir(bench):
    return ROOT / "reference" / bench / "gold"


@pytest.mark.parametrize("bench", FIELD_BENCHMARKS)
def test_field_benchmark_still_has_its_snapshot(bench):
    """The 101-step field snapshot this suite already locked at 0.0 diff."""
    d = gold_dir(bench)
    assert (d / "fault.00101.nc").is_file()
    assert (d / "fault.00101.csv").is_file()


@pytest.mark.parametrize("station", BP8_STATIONS)
def test_bp8_has_all_nine_stations(station):
    d = gold_dir("bp8")
    assert (d / f"fltst_strk{station}.csv").is_file()
    summary = json.loads((d / "summary.json").read_text())
    assert station in summary, f"summary.json carries no numbers for station {station}"


@pytest.mark.parametrize("profile", BP8_PROFILES)
def test_bp8_has_all_ten_section43_profiles(profile):
    d = gold_dir("bp8")
    assert (d / f"{profile}.csv").is_file(), (
        f"reference/bp8/gold/{profile}.csv is missing; "
        "resampleBP8Profiles.py and checkBP8Submission both exercise this "
        "file with no oracle to catch a regression in it"
    )


def test_bp8_has_a_fault_plane_snapshot():
    """BP5/BP7 have compared a fault-plane snapshot for a long time; BP8 never
    did. Filename differs from BP5/BP7 (step 05301, not 00101) because BP8 is
    a single aseismic run, not a multi-cycle one -- see reference/bp8/README.md.
    """
    d = gold_dir("bp8")
    hits = glob.glob(str(d / "fault.*.nc"))
    assert hits, "reference/bp8/gold has no fault-plane netCDF snapshot"
    csv_hits = glob.glob(str(d / "fault.*.csv"))
    assert csv_hits, "reference/bp8/gold has no flattened CSV of the fault-plane snapshot"


def test_bp8_fault_snapshot_csv_matches_the_netcdf():
    """reference/bp8/gold/fault.05301.csv is the flattened, human-eyeballable
    twin of fault.05301.nc (same convention as BP5/BP7's fault.00101.csv). A
    reader that only checks the .nc would miss the .csv going stale if either
    is regenerated without the other -- this cross-checks them row for row.
    """
    netCDF4 = pytest.importorskip("netCDF4")
    import numpy as np

    d = gold_dir("bp8")
    nc_path = d / "fault.05301.nc"
    csv_path = d / "fault.05301.csv"
    assert csv_path.is_file(), "reference/bp8/gold/fault.05301.csv is missing"

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
    d = gold_dir("bp8")
    csv_path = d / "global.csv"
    assert csv_path.is_file(), "reference/bp8/gold/global.csv is missing"

    import numpy as np
    rows = np.genfromtxt(csv_path, delimiter=",", skip_header=1)
    g = json.loads((d / "summary.json").read_text())["global"]
    assert rows[-1, 0] / 86400.0 == pytest.approx(g["t_end_d"], rel=1e-3)
    assert rows[:, 1].max() == pytest.approx(g["peak_Vmax_log10"], abs=1e-3)
