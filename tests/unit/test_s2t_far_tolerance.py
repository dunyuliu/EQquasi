"""Issue 270: the EQ_S2T far-field proportionality check (ddS2tApply) compares
constmp/w across dofs, and constmp is a running sum, so the spread is round-off
that grows with the step count n. The arithmetic of the issue is reproduced
here: the spread passes the old fixed 1e-12 yet stays below 4*eps*n, the
tolerance the solver now uses. The last test is a source-text guard (not a
run): it fails if the check goes back to a bare 1e-12."""

import re

import numpy as np

from conftest import read

EPS = np.finfo(float).eps


def _spread():
    rng = np.random.default_rng(0)
    dt = np.concatenate([6.3e5 * rng.uniform(0.5, 2, 7000), np.full(17000, 1 / 30.)])
    r = np.array([13.85e-3, 3.15e-3]) / 3.15e7
    c = np.zeros(2)
    for d in dt:
        c += r * d
    a = c / r
    return abs(a[0] - a[1]) / a.max(), dt.size


def test_spread_exceeds_old_limit_below_scaled_tolerance():
    spread, n = _spread()
    assert spread > 1.0e-12, spread            # the old check aborted here
    assert spread < 4.0 * EPS * n, (spread, 4.0 * EPS * n)
    assert max(1.0e-12, 4.0 * EPS * n) > spread


def test_source_text_guard_tolerance_scales_with_step_count():
    src = read("src/solveTimeLoopPETScDD.f90")
    assert "s2t_far_nacc" in read("src/globalvar.f90")
    assert "s2t_far_nacc = s2t_far_nacc + 1" in src
    assert "tol = max(1.0d-12, 4.0d0*epsilon(1.0d0)*dble(s2t_far_nacc))" in src
    assert not re.search(r"amax - amin > 1\.0d-12|zmax > 1\.0d-12", src)
    for f in ("src/eqquasi.f90", "src/onFaultHelpers.f90"):
        assert "s2t_far_nacc = 0" in read(f), f
