"""scripts/generateFaultInterface must build a plane that dips par.dip.

Issue 271: the per-row offset was dx*cos(dip), which builds atan(1/cos(dip))
(57.3 degrees for 50). The rows of bFault_Rough_Geometry.txt are dx apart in
depth, so the offset is dx/tan(dip). The check reads the written file back and
recovers the dip from its own y column, never from the formula under test."""

import importlib.machinery
import importlib.util
import math
import os
import sys
from types import SimpleNamespace

import numpy as np
import pytest

from conftest import ROOT

NX = NZ = 6
DX = 100.0


@pytest.fixture(scope="module")
def gen():
    os.environ.setdefault("MPLBACKEND", "Agg")
    sys.path.insert(0, str(ROOT / "scripts"))
    path = str(ROOT / "scripts" / "generateFaultInterface")
    loader = importlib.machinery.SourceFileLoader("generateFaultInterface", path)
    spec = importlib.util.spec_from_loader("generateFaultInterface", loader)
    mod = importlib.util.module_from_spec(spec)
    loader.exec_module(mod)
    return mod


def _built_dip(gen, tmp_path, dip):
    par = SimpleNamespace(fxmin=0.0, fzmin=0.0)
    gen.generateFaultInterface(par, str(tmp_path), NX, NZ, DX, DX, 1, dip, seedId=1)
    a = np.loadtxt(tmp_path / "bFault_Rough_Geometry.txt", skiprows=2)
    y = a[:, 0].reshape(NX, NZ)          # index i*nz + j: strike, then depth row
    dy = y[:, :-1] - y[:, 1:]             # y offset per dx of depth
    assert np.allclose(dy, dy[0, 0]), "planar fault must have one offset per row"
    return math.degrees(math.atan2(DX, dy[0, 0]))


@pytest.mark.parametrize("dip", [30.0, 45.0, 50.0, 60.0, 70.0, 80.0, 90.0])
def test_built_dip_equals_par_dip(gen, tmp_path, dip):
    assert _built_dip(gen, tmp_path, dip) == pytest.approx(dip, abs=1e-4), (
        f"par.dip={dip} built a plane that dips something else")


def test_dip_zero_is_refused(gen, tmp_path):
    with pytest.raises(SystemExit):
        gen.generateFaultInterface(SimpleNamespace(fxmin=0.0, fzmin=0.0),
                                   str(tmp_path), NX, NZ, DX, DX, 1, 0.0, seedId=1)
