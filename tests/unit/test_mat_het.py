"""par.mat_het helpers in scripts/defaultParameters.py: the field mirrors the
solver's shear_modulus_at (src/globalvar.f90), and a malformed dict stops
setup instead of running with a guessed value (rule 2)."""

import sys

import numpy as np
import pytest

from conftest import ROOT

sys.path.insert(0, str(ROOT / "scripts"))
from defaultParameters import (mat_het_checked, shear_modulus_at,  # noqa: E402
                               shear_modulus_range)

MU0 = 36.0e9


class Par:
    rou = 2670.0
    vs = (MU0 / 2670.0) ** 0.5


def case3b():
    p = Par()
    p.mat_het = dict(amp=10.0e9, lx=30.0e3, lz=30.0e3, z0=11.0e3,
                     box=(90.0e3, 40.0e3, 80.0e3))
    return p


def test_unset_is_uniform():
    p = Par()
    assert shear_modulus_range(p) == pytest.approx((MU0, MU0))
    assert shear_modulus_at(p, 5e3, 0.0, 7e3) == pytest.approx(MU0)


def test_field_values_and_depth_sign():
    p = case3b()
    # depth 4 km: cos(2 pi (4+11)/30) = -1; depth 19 km: cos(2 pi) = +1.
    mu = shear_modulus_at(p, [0.0, 0.0, 15e3], 0.0, [4e3, 19e3, 4e3])
    # x = 15 km flips the sign: (-1)(-1) = +1.
    assert mu == pytest.approx([26e9, 46e9, 46e9])


def test_box_cutoff():
    p = case3b()
    assert shear_modulus_at(p, 120e3, 0.0, 4e3) == pytest.approx(MU0)
    assert shear_modulus_at(p, 0.0, 50e3, 4e3) == pytest.approx(MU0)
    assert shear_modulus_at(p, 0.0, 0.0, 85e3) == pytest.approx(MU0)
    assert [x / 1e9 for x in shear_modulus_range(p)] == pytest.approx([26, 46])


@pytest.mark.parametrize("change", [
    {"Z0": 1.0},                       # typo: z0 missing, unknown key
    {"box": (90e3, 0.0, 80e3)},        # zero extent
    {"lx": 0.0},
    {"amp": 40e9},                     # |amp| >= mu0
])
def test_malformed_stops(change):
    p = case3b()
    m = dict(p.mat_het)
    if "Z0" in change:
        m.pop("z0")
    m.update(change)
    p.mat_het = m
    with pytest.raises(SystemExit):
        mat_het_checked(p)
