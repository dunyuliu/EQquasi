#! /usr/bin/env python3
"""Contour the rupture-time front over the fault plane.

Rupture time is `fnft`, written as the 16th column of `cplot_EQquasi.txt` by
`output_onfault_transfer` (src/library_output.f90:226). A node that has not
ruptured carries the sentinel -1000, not zero, so it must be masked rather
than plotted -- contouring the raw column produces a meaningless cliff at the
rupture front and a plot that looks like a result.

Columns of cplot_EQquasi.txt, in order:

    1  x            9  vym (master)
    2  -z          10  vzm
    3  slip rate   11  vxs (slave)
    4  state       12  vys
    5  tau_strike  13  vzs
    6  tau_dip     14  (unused)
    7  sigma_eff   15  (unused)
    8  vxm         16  fnft, rupture time (s), -1000 = never ruptured

Applicability, which is worth knowing before you go looking for a front:

  BP5, BP5-dip90   seismic; a full cycle produces a real front. The 101-step
                   test compsets only reach artificial nucleation, so expect a
                   few tens of ruptured nodes, not a propagating front.
  BP7              seismic, but nucleates later than the test compset runs.
  BP8              aseismic by construction -- velocity-strengthening
                   everywhere, a > b, so nothing ever ruptures and fnft stays
                   at the sentinel. An empty plot here is the correct answer,
                   not a failure.

Multi-fault runs: cplot_EQquasi.txt stacks every fault's nodes in one file;
nodes are regridded from their coordinates, so all faults appear on the one
map (they occupy disjoint strike ranges).

Reads:  cplot_EQquasi.txt
Writes: rupture_time.png per cycle
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import plotutils as pu

pu.apply_style()
import matplotlib.pyplot as plt
import numpy as np

PATTERNS = ["cplot_EQquasi.txt"]
UNRUPTURED = -1000.0
FNFT_COL = 15          # zero-based


def process_dir(rdir, outdir, levels):
    a = np.loadtxt(os.path.join(rdir, "cplot_EQquasi.txt"))
    x, z, fnft = a[:, 0], a[:, 1], a[:, FNFT_COL]

    ruptured = fnft > UNRUPTURED / 2.0        # well clear of the sentinel
    n = int(ruptured.sum())
    print(f"  {n} of {len(fnft)} fault nodes ruptured "
          f"({100.0 * n / len(fnft):.1f} %)")
    if n == 0:
        print("  nothing ruptured -- for BP8 that is expected (aseismic by "
              "construction); for BP5/BP7 the run has not nucleated yet.")
        return
    print(f"  rupture time: {fnft[ruptured].min():.4g} to "
          f"{fnft[ruptured].max():.4g} s")

    # Reshape onto the fault grid. Nodes are written in a single loop, so infer
    # the grid from the unique coordinates rather than assuming an ordering.
    xs, zs = np.unique(x), np.unique(z)
    grid = np.full((len(zs), len(xs)), np.nan)
    grid[np.searchsorted(zs, z), np.searchsorted(xs, x)] = \
        np.where(ruptured, fnft, np.nan)

    fig, ax = plt.subplots(figsize=(10, 5.5))
    m = ax.contourf(xs / 1e3, zs / 1e3, grid, levels=levels, cmap="turbo")
    cs = ax.contour(xs / 1e3, zs / 1e3, grid, levels=levels,
                    colors="k", linewidths=0.4, alpha=0.5)
    ax.clabel(cs, inline=True, fontsize=7, fmt="%.2g")
    cb = plt.colorbar(m, ax=ax)
    cb.set_label("rupture time (s)")
    ax.set_xlabel("along strike (km)")
    ax.set_ylabel("z (km)")
    ax.set_aspect("equal")
    ax.set_title(f"Rupture-time contours — "
                 f"{os.path.basename(os.path.abspath(rdir))}\n"
                 f"{n}/{len(fnft)} nodes ruptured; unruptured masked")
    fig.tight_layout()
    pu.save(fig, pu.out_path(rdir, "rupture_time.png", outdir), dpi=150)


def main():
    ap = pu.make_parser(__doc__, "plotRuptureTime.py",
                        "rupture_time.png into each cycle's directory "
                        "(empty run: a message, no figure).", PATTERNS)
    ap.add_argument("--levels", type=int, default=20,
                    help="number of contour levels (default 20)")
    args = ap.parse_args()
    for label, rdir in pu.resolve_targets(args.dirs, PATTERNS,
                                          "plotRuptureTime.py"):
        print(f"processing {label or rdir} ({rdir})")
        process_dir(rdir, args.outdir, args.levels)
    return 0


if __name__ == "__main__":
    sys.exit(main())
