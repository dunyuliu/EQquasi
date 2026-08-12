"""Guard: every gold fault snapshot must have the shape the writer produces.

`netcdf_write_on_fault` (src/netcdf_io.f90) carries an explicit `nid_fault`
dimension of size `ntotft`, so on-fault variables are stored 3-D as
(nid_fault, nid_dip, nid_strike). Gold predating that change is 2-D.

The `mfault` merge introduced the dimension and left the BP5, BP5-dip90 and BP7
gold at the old 2-D shape. CI caught it -- `check.test.py` compares with
xarray's `identical()`, which rejects differing dimensions -- but
`scripts/plotAgainstGold.py` did not: its `load_field` keys on
(variable, fault index) and maps (nz, nx) and (1, nz, nx) onto the same key, so
the shape difference was consumed before any comparison and it reported
"matches gold" for fourteen hours across four releases.

This test is the cheap version of that CI signal: it needs no build and no run,
so a stale gold shape fails in seconds rather than after an eleven-minute
pipeline. It checks the stored dimensions only -- values are the e2e tier's job.
"""

import glob
import os

import pytest

from conftest import ROOT

pytestmark = pytest.mark.contract

FAULT_SNAPSHOTS = sorted(glob.glob(str(ROOT / "reference" / "*" / "fault.*.nc")))

# Variables are (nid_fault, nid_dip, nid_strike). Coordinate variables named
# nid_* are 1-D by construction and are not on-fault fields.
REQUIRED_DIMS = {"nid_fault", "nid_dip", "nid_strike"}


@pytest.mark.skipif(not FAULT_SNAPSHOTS, reason="no gold fault snapshots")
@pytest.mark.parametrize("path", FAULT_SNAPSHOTS,
                         ids=[os.path.relpath(p, str(ROOT)) for p in FAULT_SNAPSHOTS])
def test_gold_snapshot_carries_the_fault_dimension(path):
    nc = pytest.importorskip("netCDF4")
    d = nc.Dataset(path)
    dims = set(d.dimensions)
    missing = REQUIRED_DIMS - dims
    assert not missing, (
        f"{os.path.relpath(path, str(ROOT))} is missing {sorted(missing)}; it has "
        f"{sorted(dims)}. netcdf_write_on_fault always writes nid_fault, so this "
        "gold predates that change and no current run can match it. Regenerate it."
    )


@pytest.mark.skipif(not FAULT_SNAPSHOTS, reason="no gold fault snapshots")
@pytest.mark.parametrize("path", FAULT_SNAPSHOTS,
                         ids=[os.path.relpath(p, str(ROOT)) for p in FAULT_SNAPSHOTS])
def test_gold_snapshot_variables_are_three_dimensional(path):
    nc = pytest.importorskip("netCDF4")
    d = nc.Dataset(path)
    flat = [v for v, var in d.variables.items()
            if not v.startswith("nid_") and var.ndim != 3]
    assert not flat, (
        f"{os.path.relpath(path, str(ROOT))}: {flat} are not "
        "(nid_fault, nid_dip, nid_strike). A 2-D variable here compares clean "
        "against a 3-D run under any comparator that reduces per fault, which is "
        "how this went unnoticed once already."
    )
