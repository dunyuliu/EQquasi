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


def test_fortran_routine_itself_matches_the_mirror(tmp_path):
    """Compile the REAL subroutine out of src/fric.f90 and compare numbers, so a flipped branch, a swapped
    argument or a dropped guard in the Fortran fails here (the string checks above cannot see those)."""
    import shutil, subprocess, pytest
    fc = shutil.which("gfortran")
    if fc is None:
        pytest.skip("no gfortran")
    src = read("src/fric.f90")
    body = src[src.index("SUBROUTINE lapusta_xi"):src.index("end SUBROUTINE lapusta_xi") + len("end SUBROUTINE lapusta_xi")]
    cases = [(kL(36e9, 1000.0), A, B, 0.5), (kL(46e9, 1000.0), A, B, 0.5), (kL(46e9, 500.0), A, B, 0.5),
             (kL(36e9, 1000.0), A, B, 0.05), (1.2 * (B - A), A, B, 0.5), (3.0 * A, A, 0.5 * A, 0.5),
             (kL(36e9, 1000.0), 0.0, B, 0.5), (0.5 * (B - A), A, B, 0.5), (kL(36e9, 1000.0), B, A, 0.5)]
    drv = "module globalvar\n integer, parameter :: dp = kind(1.0d0)\nend module\n" + body + "\nprogram t\n use globalvar, only : dp\n real(dp) :: xi\n"
    for c in cases:
        drv += " call lapusta_xi(%s, xi); print '(es24.16)', xi\n" % ", ".join("%.17gd0" % v for v in c)
    drv += "end program\n"
    (tmp_path / "t.f90").write_text(drv)
    subprocess.run([fc, "-ffree-line-length-none", "-o", "t", "t.f90"], cwd=tmp_path, check=True)
    out = [float(x) for x in subprocess.run(["./t"], cwd=tmp_path, check=True, capture_output=True, text=True).stdout.split()]
    assert len(out) == len(cases)
    for c, f in zip(cases, out):
        m = lapusta_xi(*c)
        assert (f <= 0 and m <= 0) or f == pytest.approx(m, rel=1e-12), (c, f, m)


def test_call_site_passes_the_right_quantities():
    f = " ".join(read("src/faulting.f90").split())
    assert ("call lapusta_xi(xi_gamma*shear_modulus_at(x(1,isn), x(2,isn), x(3,isn))/dx*fric(FR_RSF_DC,i,ift), & "
            "fric(FR_RSF_A,i,ift)*abs(fric(FR_TNRM,i,ift)), fric(FR_RSF_B,i,ift)*abs(fric(FR_TNRM,i,ift)), & ksi, xi_node)") in f


def test_audit_fixes_pinned():
    """PR #72 audit: NaN xi_gamma refused (Python and Fortran), per-node xi only at
    frictional nodes, and a failed per-node rule aborts every rank."""
    f = read("src/faulting.f90")
    assert "call MPI_Abort(MPI_COMM_WORLD, 13, ierr_abort)" in f
    assert "stop 13" not in f
    assert ("if (xi_gamma > 0.0d0 .and. abs(x(3,isn)) <= abs(zminc) &" in f
            and ".and. x(1,isn) <= xmaxc .and. x(1,isn) >= xminc) then" in f)
    assert ".not. (xi_gamma >= 0.0d0 .and. xi_gamma <= 1.0d6)" in read("src/read_input.f90")
    assert "math.isfinite(float(par.xi_gamma))" in read("scripts/case.setup")
