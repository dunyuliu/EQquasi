"""The fluid domain must be Omega_f = [-l_f, l_f]^2, not (N*dx)^2.

BP8 confines fluid flow to the frictional domain and imposes zero flux across
its boundary (eq. 18). Omega_f is the closed square [-400, 400]^2, so a node
grid with N points per side must store fluid in ((N-1)*dx)^2 = 640000 m^2.

Through v1.4.6 every node carried a full dx^2 cell, including boundary nodes, so
the storage area was (N*dx)^2 -- 722500 m^2 at dx = 50 m, 12.9 % too much. The
injected fluid was spread through too much rock and the late-time uniform
pressure came out 11.4 % low: 1.4948 MPa measured against 1.6875 MPa required.
The reference reads ~1.7 MPa.

The mechanism is the discrete Laplacian. Dropping a neighbour that lies outside
Omega_f, without rescaling, is not a zero-flux condition -- it models a full cell
hanging off the edge. The zero-flux condition is a mirrored ghost node, which
doubles the normal-direction term; equivalently, boundary nodes carry half a
cell (a quarter at a corner) and the flux divergence is divided by that fraction.

Two guards: the scheme itself is checked numerically against the exact
conservation answer, and the Fortran is checked for the terms that implement it.
"""

import math
import re

import numpy as np
import pytest

from conftest import read, strip_fortran_comments

# Table 1
DX_CASES = [(50.0, 17), (10.0, 81)]
L_F = 400.0
ALPHA, BETA, PHI, LFWID = 0.05, 1.0e-8, 0.1, 1.0
Q0, TOFF, LGAUSS = 0.003, 100 * 3600.0, 50.0


def cell_fractions(n):
    fx = np.ones((n, n))
    fz = np.ones((n, n))
    fx[0, :] = fx[-1, :] = 0.5
    fz[:, 0] = fz[:, -1] = 0.5
    return fx, fz


@pytest.mark.parametrize("dx,n", DX_CASES)
def test_storage_area_equals_omega_f(dx, n):
    """sum of cell areas must be exactly (2*l_f)^2, independent of dx."""
    fx, fz = cell_fractions(n)
    area = (dx * dx * fx * fz).sum()
    assert area == pytest.approx((2 * L_F) ** 2, rel=1e-12), (
        f"dx={dx}: storage area {area:.1f} m^2, required {(2*L_F)**2:.1f} m^2. "
        f"Full cells everywhere would give {(n*dx)**2:.1f}."
    )


@pytest.mark.parametrize("dx,n", DX_CASES)
def test_late_time_pressure_matches_exact_conservation(dx, n):
    """After shutoff the zero-flux boundary traps everything injected.

    p_uniform = Q0*t_off / (beta*phi*A*L_fwid), with A the true Omega_f area.
    This exercises the operator, the source normalisation and the boundary
    condition together, and is independent of the Green's function.
    """
    xs = np.linspace(-L_F, L_F, n)
    X, Z = np.meshgrid(xs, xs, indexing="ij")
    fx, fz = cell_fractions(n)
    a = dx * dx * fx * fz

    w = np.exp(-(X**2 + Z**2) / (2.0 * LGAUSS**2))
    w /= w.sum()

    p = np.zeros((n, n))
    dt = 0.2 * dx * dx / ALPHA
    nst = int(math.ceil(TOFF / dt))
    dt = TOFF / nst

    def lap(pp):
        lx = np.zeros_like(pp)
        lz = np.zeros_like(pp)
        lx[1:, :] += pp[:-1, :] - pp[1:, :]
        lx[:-1, :] += pp[1:, :] - pp[:-1, :]
        lz[:, 1:] += pp[:, :-1] - pp[:, 1:]
        lz[:, :-1] += pp[:, 1:] - pp[:, :-1]
        return (lx / fx + lz / fz) / (dx * dx)

    src = Q0 * w / (BETA * PHI * a * LFWID)
    for _ in range(nst):
        p = p + dt * (ALPHA * lap(p) + src)
    # shut off, then run long enough to homogenise
    for _ in range(40 * nst):
        p = p + dt * ALPHA * lap(p)

    expected = Q0 * TOFF / (BETA * PHI * (2 * L_F) ** 2 * LFWID)
    assert p.mean() == pytest.approx(expected, rel=1e-6), (
        f"dx={dx}: late-time uniform pressure {p.mean()/1e6:.4f} MPa, "
        f"expected {expected/1e6:.4f} MPa"
    )
    assert (p.max() - p.min()) < 1.0, "field failed to homogenise"


def test_fortran_scales_the_laplacian_by_the_cell_fraction():
    """The Fortran must divide each direction by its own cell fraction."""
    src = strip_fortran_comments(read("src/porepressure.f90"))
    m = re.search(r"lap\s*=\s*\(\s*lapx\s*/\s*pfFx\(i\)\s*\+\s*lapz\s*/\s*pfFz\(i\)\s*\)"
                  r"\s*/\s*\(\s*dx\s*\*\s*dx\s*\)", src)
    assert m, (
        "porepressure.f90 no longer divides the Laplacian by pfFx/pfFz. Without "
        "that, dropping an out-of-domain neighbour models a full cell outside "
        "Omega_f instead of a zero-flux boundary, and Omega_f becomes (N*dx)^2."
    )


def test_fortran_source_uses_the_node_cell_area():
    """sum_i src_i * A_i must be q_inj/(beta*phi*L_fwid), so src divides by A_i."""
    src = strip_fortran_comments(read("src/porepressure.f90"))
    assert re.search(r"acell\s*=\s*area\s*\*\s*pfFx\(i\)\s*\*\s*pfFz\(i\)", src), \
        "the per-node cell area acell is no longer computed"
    assert "fluid_phi * area * fluid_Lfwid" not in src, (
        "the source still divides by the full dx^2 cell rather than acell; "
        "the injected volume will not be conserved on boundary cells"
    )
