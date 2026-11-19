"""Guard: every gold fault snapshot must have the shape the writer produces.

`netcdf_write_on_fault` (src/netcdf_io.f90) carries an explicit `nid_fault`
dimension of size `ntotft`, so on-fault variables are stored 3-D as
(nid_fault, nid_dip, nid_strike). Gold predating that change is 2-D.

The `mfault` merge introduced the dimension and left the BP5, BP5-dip90 and BP7
gold at the old 2-D shape. CI caught it -- `tests/e2e/test_benchmarks.py` compares with
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

FAULT_SNAPSHOTS = sorted(glob.glob(str(ROOT / "data" / "*" / "fault.*.nc")))

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


# ---------------------------------------------------------------------------
# Storage width. The NetCDF writers store 8-byte reals since 1.24.8 (src/netcdf_io.f90
# follows dp; owner decision, board row 1, 2026-11-18), and the tier bar for those
# files is 1e-10. A 4-byte file under data/ in a regenerated set would make that bar
# unreachable (float32 holds ~6e-8), and nothing else checked the width of the
# reference files. A set joins REGENERATED_8BYTE in the commit that replaces its
# files, never before: a set still holding 4-byte gold is the next regeneration's
# job, not a failure of this one.
# ---------------------------------------------------------------------------

REGENERATED_8BYTE = (
    "bp1002.qdc.2500/cycle0",
    "test.bp5.qdc.2000/cycle0-step101-fast",
    "test.bp5.qdc.dip90.2000/cycle0",
    "test.stepover.qdc.1000/cycle0",
    "test.stepover.qdc.con.1000/cycle0",
    "test.bp7.qdc.a.10/cycle0-step101-fast",
    "test.bp7.qdc.a.10/cycle0",
    "bp5.qdc.kink.2000/cycle0",
    "test.bp8.qdc.gs.10",
    "test.bp5.qdc.2000/cycle0",
)
REGENERATED_FILES = sorted(p for s in REGENERATED_8BYTE
                           for p in glob.glob(str(ROOT / "data" / s / "*.nc")))


def assert_eight_byte_reals(path):
    """Raise AssertionError naming every real variable of `path` not stored in 8 bytes.

    A file with no real variable at all also raises: a width check that passes
    because nothing was checked is the failure this guard exists to prevent.
    """
    nc = pytest.importorskip("netCDF4")
    d = nc.Dataset(path)
    real = {v: var.dtype.itemsize for v, var in d.variables.items() if var.dtype.kind == "f"}
    d.close()
    assert real, f"{path} holds no real variable; the width check checked nothing"
    narrow = {v: w for v, w in real.items() if w != 8}
    assert not narrow, (
        f"{path}: {narrow} are stored as {sorted(set(narrow.values()))}-byte reals (4-byte means float32); "
        "the regenerated sets hold 8-byte reals (1.24.8). Regenerate the file with the 1.24.8 binary, "
        "or take the set out of REGENERATED_8BYTE in the commit that decides to keep it."
    )


@pytest.mark.parametrize("subdir", REGENERATED_8BYTE)
def test_regenerated_set_is_present(subdir):
    found = glob.glob(str(ROOT / "data" / subdir / "fault.*.nc"))
    assert found, f"data/{subdir} has no fault.*.nc; REGENERATED_8BYTE names a set that is not there"


@pytest.mark.parametrize("path", REGENERATED_FILES,
                         ids=[os.path.relpath(p, str(ROOT)) for p in REGENERATED_FILES])
def test_regenerated_gold_stores_8_byte_reals(path):
    assert_eight_byte_reals(path)


def test_width_check_fails_loudly_on_a_4_byte_file(tmp_path):
    nc = pytest.importorskip("netCDF4")
    for code, ok in (("f4", False), ("f8", True)):
        p = str(tmp_path / f"w_{code}.nc")
        d = nc.Dataset(p, "w")
        d.createDimension("n", 3)
        v = d.createVariable("slip", code, ("n",))
        v[:] = [1.0, 2.0, 3.0]
        d.close()
        if ok:
            assert_eight_byte_reals(p)
        else:
            with pytest.raises(AssertionError, match="4-byte"):
                assert_eight_byte_reals(p)


def test_width_check_fails_on_a_file_with_no_real_variable(tmp_path):
    nc = pytest.importorskip("netCDF4")
    p = str(tmp_path / "ints.nc")
    d = nc.Dataset(p, "w")
    d.createDimension("n", 3)
    d.createVariable("nid", "i4", ("n",))[:] = [1, 2, 3]
    d.close()
    with pytest.raises(AssertionError, match="checked nothing"):
        assert_eight_byte_reals(p)
