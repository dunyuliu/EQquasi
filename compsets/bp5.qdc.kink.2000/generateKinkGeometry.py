#! /usr/bin/env python3
"""Write bFault_Rough_Geometry.txt for a strike-slip fault with a kink.

Liu, Duan & Luo (2020, GJI 220, 598-609) model a 60 km strike-slip fault with
a 10 degree bend at its centre. EQquasi has no bent-fault geometry: faultgeom
places each fault on a constant-y plane, so a segment oblique to x cannot be
expressed that way. The rough-fault path can express it. insert_rough_fault
(src/func_lib.f90) displaces a node's y by a value read per (x, z) from
bFault_Rough_Geometry.txt, so any single-valued y(x, z) is available -- a
fractal surface in the rough case, and here a piecewise-linear kink.

The file is the same one script/generateFaultInterface writes, and the
format is fixed by read_fault_rough_geometry (src/read_input.f90):

    line 1        nx   nz   0
    line 2        dx   fxmin   fzmin
    lines 3..     y    dy/dx   dy/dz        x outer, z inner

y is the fault-normal offset of the surface; the two gradients are read into
pfx and pfz. This writes all three, with the gradients taken from the same
array by np.gradient so they stay consistent with y rather than being derived
analytically and drifting from it at the kink.

Geometry, following the paper's figure 1(b): the left segment lies along the
x axis, and the right segment leaves the bend at `angle` degrees, so

    y(x) = 0                                for x <= xBend
    y(x) = (x - xBend) * tan(angle)         for x >  xBend

The kink is vertical: y does not vary with z, so both segments are vertical
planes and the bend is purely in map view, as in the paper.

Usage
-----
    ./generateKinkGeometry.py [case_dir]

Reads par from the case's user_defined_params.py (fxmin, fzmin, dx, nfx,
nfz, and the kink parameters below), writes bFault_Rough_Geometry.txt and a
PNG of the resulting trace beside it. case.setup calls
script/generateFaultInterface, not this script, so a compset that wants a
kink sets par.insertFaultType = 0 and runs this instead -- see the compset's
own comment.

Kink parameters are module-level constants in the compset, not attributes on
par -- par is the solver's schema and the solver never reads these:
    KINK_ANGLE_DEG   bend angle in degrees
    KINK_X           along-strike position, m
"""

import os
import sys
from math import pi, tan

import numpy as np


def kink_surface(x, z, kink_x, angle_deg):
    """y(x, z) for a vertical fault kinked by `angle_deg` at x = kink_x.

    Returned shape is (nz, nx), matching the (j, i) indexing the file format
    and insert_rough_fault both use.
    """
    y = np.where(x > kink_x, (x - kink_x) * tan(angle_deg / 180.0 * pi), 0.0)
    return np.tile(y, (len(z), 1))


def write_geometry(case_dir, fxmin, fzmin, dx, dz, nx, nz,
                   kink_x=0.0, angle_deg=10.0):
    """Write bFault_Rough_Geometry.txt; return the surface for inspection."""
    x = fxmin + np.arange(nx) * dx
    z = fzmin + np.arange(nz) * dz
    surface = kink_surface(x, z, kink_x, angle_deg)

    # Gradients are taken from the surface itself, not from tan(angle): at the
    # kink the analytic derivative is undefined, and the solver needs the same
    # discrete slope the mesher will see.
    dsdx = np.gradient(surface / dx, axis=1)
    dsdz = np.gradient(surface / dz, axis=0)

    if np.abs(dsdx).max() > 0.2 or np.abs(dsdz).max() > 0.2:
        print(f"  warning: max|dy/dx| = {np.abs(dsdx).max():.3f}, "
              f"max|dy/dz| = {np.abs(dsdz).max():.3f}; over 0.2 the element "
              "distortion starts to matter (generateFaultInterface warns at "
              "the same threshold)")

    out = np.zeros((nx * nz + 2, 3))
    out[0, 0], out[0, 1] = nx, nz
    out[1, 0], out[1, 1], out[1, 2] = dx, fxmin, fzmin
    n = 2
    for i in range(nx):
        for j in range(nz):
            out[n, 0] = surface[j, i]
            out[n, 1] = dsdx[j, i]
            out[n, 2] = dsdz[j, i]
            n += 1

    path = os.path.join(case_dir, "bFault_Rough_Geometry.txt")
    np.savetxt(path, out, fmt="%f", delimiter="\t")
    print(f"  wrote {path}")
    print(f"  kink {angle_deg} deg at x = {kink_x/1e3:g} km; "
          f"y spans {surface.min():.1f} .. {surface.max():.1f} m over "
          f"x = {x[0]/1e3:g} .. {x[-1]/1e3:g} km")
    return x, z, surface


def plot(case_dir, x, z, surface):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("  matplotlib not available; skipping the figure")
        return
    fig, ax = plt.subplots(1, 2, figsize=(13, 4))
    ax[0].plot(x / 1e3, surface[-1] / 1e3, lw=2)
    ax[0].set_xlabel("along strike, x (km)")
    ax[0].set_ylabel("fault-normal, y (km)")
    ax[0].set_title("Fault trace in map view")
    ax[0].set_aspect("equal")
    ax[0].grid(alpha=0.3)
    m = ax[1].pcolormesh(x / 1e3, z / 1e3, surface / 1e3, shading="nearest")
    plt.colorbar(m, ax=ax[1], label="y (km)")
    ax[1].set_xlabel("along strike, x (km)")
    ax[1].set_ylabel("z (km)")
    ax[1].set_title("y(x, z) written to the interface file")
    fig.tight_layout()
    out = os.path.join(case_dir, "cKinkGeometry.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"  wrote {out}")


def main():
    # The argument is where the geometry is WRITTEN -- input/ in a case laid
    # out by create.newcase. user_defined_params.py is not there: it lives in
    # the case root, one level up, and tool/ holds what it imports. Python
    # puts this script's own directory on sys.path, not the caller's, so
    # running `python3 input/generateKinkGeometry.py input` from the case root
    # failed with ModuleNotFoundError until both were added explicitly.
    out_dir = sys.argv[1] if len(sys.argv) > 1 else "."
    case_dir = sys.argv[2] if len(sys.argv) > 2 else (
        os.path.dirname(os.path.abspath(out_dir)) or ".")
    for d in (case_dir, os.path.join(case_dir, "tool"), os.path.abspath(out_dir)):
        if os.path.isdir(d):
            sys.path.insert(0, os.path.abspath(d))
    import user_defined_params as udp
    par = udp.par
    case_dir = out_dir

    nx = int(round((par.fxmax - par.fxmin) / par.dx)) + 1
    nz = int(round((par.fzmax - par.fzmin) / par.dz)) + 1
    print(f"generateKinkGeometry: {nx} x {nz} points at dx = {par.dx:g} m")
    x, z, s = write_geometry(
        case_dir, par.fxmin, par.fzmin, par.dx, par.dz, nx, nz,
        kink_x=udp.KINK_X, angle_deg=udp.KINK_ANGLE_DEG)
    plot(case_dir, x, z, s)
    return 0


if __name__ == "__main__":
    sys.exit(main())
