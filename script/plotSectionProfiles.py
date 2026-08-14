#! /usr/bin/env python3
"""Plot the section 4.3 evolution profiles that go into a BP8 submission.

These ten ASCII files are two thirds of a BP8 upload and nothing plotted them,
so a result could be submitted without anyone having looked at most of it:

    slip_2_strike.dat          slip_2_depth.dat
    slip_3_strike.dat          slip_3_depth.dat
    shear_stress_2_strike.dat  shear_stress_2_depth.dat
    shear_stress_3_strike.dat  shear_stress_3_depth.dat
    pore_pressure_strike.dat   pore_pressure_depth.dat

Not the same thing as `plotProfile.py` or `plotProfiles`, which read fault-plane
netCDF snapshots. These are the section 4.3 ASCII format:

    row 0        0, 0, then the node coordinates in metres
    rows 1..N    time (s), max slip rate, then the field at each node

so each row is one profile in space and the file is its evolution in time.
Drawn as 2-D maps of distance against time -- one figure per quantity per
section, ten in all -- which is how the benchmark's own comparison pages present
them and which shows when things happen, not just what the profiles look like.

Works on a run directory (native node spacing, e.g. 17 nodes at dx = 50 m) or on
a resampled submission directory (81 nodes at 10 m) -- the node count is read
from row 0, never assumed.

Usage:
    cd <case> && plotSectionProfiles.py          all cycles
    plotSectionProfiles.py <dir>                 one directory
    plotSectionProfiles.py <dir> --lines 40      how many times to draw
"""

import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import plotutils as pu
from seasio import read_array

PATTERNS = ("*_strike.dat",)

# (stem, axis label). Ten files: five quantities on two section lines.
QUANTITIES = [
    ("slip_2", "slip$_2$ (m)"),
    ("slip_3", "slip$_3$ (m)"),
    ("shear_stress_2", r"$\tau_2$ (MPa)"),
    ("shear_stress_3", r"$\tau_3$ (MPa)"),
    ("pore_pressure", "pore pressure (MPa)"),
]
LINES = [("strike", "along strike, $x_2$ (m)"),
         ("depth", "along depth, $x_3$ (m)")]


def load_profile(path):
    """(coords, times, values) from a section 4.3 file.

    Row 0 carries the node coordinates after two leading zeros; every later row
    is time, max slip rate, then one value per node.
    """
    a = read_array(path)
    if a.ndim != 2 or a.shape[0] < 2:
        return None
    return a[0, 2:], a[1:, 0], a[1:, 2:]


def plot(rdir, outdir, nlines):
    """One figure per quantity per section: a 2-D map of distance against time.

    A family of overlaid curves shows the shape of each profile but hides when
    things happen; with 663 time slices it is also a solid block. Drawn as a
    map with time up the vertical axis, the whole evolution is one panel and
    features -- the pressure front spreading, injection switching off at
    t_off -- are read directly off the time axis.
    """
    made = 0
    for ln, xlab in LINES:
        for q, lab in QUANTITIES:
            path = os.path.join(rdir, f"{q}_{ln}.dat")
            if not os.path.exists(path):
                continue
            got = load_profile(path)
            if got is None:
                print(f"  {q}_{ln}.dat: no data rows")
                continue
            coords, times, vals = got
            td = times / 86400.0

            # Signed fields get a diverging map centred on zero, so the sign is
            # readable; strictly positive ones get a sequential map.
            lo, hi = float(np.nanmin(vals)), float(np.nanmax(vals))
            if lo < 0.0 < hi:
                m = max(abs(lo), abs(hi))
                cmap, vmin, vmax = "bwr", -m, m
            else:
                cmap, vmin, vmax = "magma", lo, hi

            fig, ax = plt.subplots(figsize=(8, 6), constrained_layout=True)
            im = ax.pcolormesh(coords, td, vals, shading="nearest",
                               cmap=cmap, vmin=vmin, vmax=vmax)
            cb = fig.colorbar(im, ax=ax, pad=0.02)
            cb.set_label(lab)
            ax.set_xlabel(xlab)
            ax.set_ylabel("time (days)")
            ax.set_title(f"{q.replace('_', ' ')}, {ln} section\n"
                         f"{os.path.basename(os.path.abspath(rdir))} — "
                         f"{len(coords)} nodes, {len(times)} times")
            pu.save(fig, pu.out_path(rdir, f"{q}_{ln}_xt.png", outdir))
            made += 1
    if not made:
        print(f"  no section 4.3 profile files in {rdir}")


def main():
    ap = pu.make_parser(__doc__, "plotSectionProfiles.py",
                        "one <quantity>_<section>_xt.png per quantity per section "
                        "(10 figures) into each results directory.", PATTERNS)
    ap.add_argument("--lines", type=int, default=30, metavar="N",
                    help=argparse.SUPPRESS)
    args = ap.parse_args()
    for label, rdir in pu.resolve_targets(args.dirs, PATTERNS,
                                          "plotSectionProfiles.py"):
        print(f"processing {label} ({rdir})")
        plot(rdir, args.outdir, args.lines)
    return 0


if __name__ == "__main__":
    sys.exit(main())
