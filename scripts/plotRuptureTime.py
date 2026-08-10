#! /usr/bin/env python3
"""Contour the rupture-time front over the fault plane.

Rupture time is `fnft`, written as the 16th column of `cplot_EQquasi.txt` by
`output_onfault_transfer` (src/library_output.f90:226). A node that has not
ruptured carries the sentinel **-1000**, not zero, so it must be masked rather
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
  BP8              **aseismic by construction** -- velocity-strengthening
                   everywhere, a > b, so nothing ever ruptures and fnft stays at
                   the sentinel. An empty plot here is the correct answer, not a
                   failure.

Usage:
    plotRuptureTime.py <run_dir> [-o out.png] [--levels N]
"""

import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

UNRUPTURED = -1000.0
FNFT_COL = 15          # zero-based


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dir")
    ap.add_argument("-o", "--out", default=None)
    ap.add_argument("--levels", type=int, default=20)
    args = ap.parse_args()

    path = os.path.join(args.run_dir, "cplot_EQquasi.txt")
    if not os.path.exists(path):
        print(f"no cplot_EQquasi.txt in {args.run_dir}")
        return 2
    a = np.loadtxt(path)
    x, z, fnft = a[:, 0], a[:, 1], a[:, FNFT_COL]

    ruptured = fnft > UNRUPTURED / 2.0        # well clear of the sentinel
    n = int(ruptured.sum())
    print(f"{n} of {len(fnft)} fault nodes ruptured "
          f"({100.0*n/len(fnft):.1f} %)")
    if n == 0:
        print("nothing ruptured -- for BP8 that is expected (aseismic by "
              "construction); for BP5/BP7 the run has not nucleated yet.")
        return 0
    print(f"rupture time: {fnft[ruptured].min():.4g} to {fnft[ruptured].max():.4g} s")

    # Reshape onto the fault grid. Nodes are written in a single loop, so infer
    # the grid from the unique coordinates rather than assuming an ordering.
    xs, zs = np.unique(x), np.unique(z)
    grid = np.full((len(zs), len(xs)), np.nan)
    ix = np.searchsorted(xs, x)
    iz = np.searchsorted(zs, z)
    grid[iz, ix] = np.where(ruptured, fnft, np.nan)

    fig, ax = plt.subplots(figsize=(9, 5))
    m = ax.contourf(xs / 1e3, zs / 1e3, grid, levels=args.levels, cmap="turbo")
    cs = ax.contour(xs / 1e3, zs / 1e3, grid, levels=args.levels,
                    colors="k", linewidths=0.4, alpha=0.5)
    ax.clabel(cs, inline=True, fontsize=6, fmt="%.2g")
    cb = plt.colorbar(m, ax=ax)
    cb.set_label("rupture time (s)")
    ax.set_xlabel("along strike (km)")
    ax.set_ylabel("depth (km)")
    ax.set_title(f"Rupture-time contours — {os.path.basename(os.path.normpath(args.run_dir))}\n"
                 f"{n}/{len(fnft)} nodes ruptured; unruptured masked", fontsize=10)
    fig.tight_layout()
    out = args.out or os.path.join(args.run_dir, "rupture_time.png")
    fig.savefig(out, dpi=130)
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
