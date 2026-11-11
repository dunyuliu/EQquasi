"""par.xi_gamma (v1.23.0): the per-node time-step factor of Lapusta et al.
(2000), eqs. 15a-c. The Python mirror is checked against hand values, and the
Fortran routine is checked to be the same expressions, so the two cannot
drift; case.setup must write the keyword line read_input.f90 looks for."""

import math
import sys

import pytest

from conftest import ROOT, read

sys.path.insert(0, str(ROOT / "scripts"))
from defaultParameters import lapusta_xi  # noqa: E402

A, B, DC = 0.005 * 25e6, 0.010 * 25e6, 0.05      # preprint Case 3 friction
G = math.pi / 4


def kL(mu, dx):
    return G * mu / dx * DC


def test_hand_values_chi_positive():
    # chi > 0: xi = A / (kL - (B - A)); stiffer rock and smaller cells -> smaller xi
    assert lapusta_xi(kL(36e9, 1000.0), A, B, 0.5) == pytest.approx(0.0970, abs=2e-4)
    assert lapusta_xi(kL(46e9, 1000.0), A, B, 0.5) == pytest.approx(0.0743, abs=2e-4)
    assert lapusta_xi(kL(46e9, 500.0), A, B, 0.5) == pytest.approx(0.0358, abs=2e-4)


def test_cap_applies():
    assert lapusta_xi(kL(46e9, 1000.0), A, B, 0.05) == 0.05


def test_chi_negative_branch():
    # kL/A = 2, (B-A)/A = 1: chi = 1/4 - 2 < 0, xi = 1 - (B-A)/kL = 1/2
    assert lapusta_xi(2 * A, A, 2 * A, 1.0) == pytest.approx(0.5)


def test_no_stable_step_is_reported_not_replaced():
    assert lapusta_xi(kL(36e9, 1000.0), 0.0, B, 0.5) <= 0       # no direct effect
    assert lapusta_xi(0.5 * (B - A), A, B, 0.5) <= 0            # cell larger than h*


def test_fortran_mirror_and_keyword_line():
    f = read("src/fric.f90")
    for expr in ("chi = 0.25d0*(kL/A - (B - A)/A)**2 - kL/A",
                 "xi = min(A/(kL - (B - A)), xic)",
                 "xi = min(1.0d0 - (B - A)/kL, xic)"):
        assert expr in f, expr
    assert 'f.write("xi_gamma  " + str(par.xi_gamma)' in read("scripts/case.setup")
    assert "index(adjustl(capsLine), 'xi_gamma') == 1" in read("src/read_input.f90")
