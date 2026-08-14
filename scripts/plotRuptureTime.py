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

Multi-fault runs: --fault selects the fault. Each fault has its OWN file --
fault 1 writes cplot_EQquasi.txt, later faults cplot_ft<N>_EQquasi.txt -- so
one figure shows one fault. This docstring used to claim the file stacked
every fault's nodes together; it never did. Before v1.13.0 the solver wrote
fault 1 only, and a cycle whose event was on another fault reported "nothing
ruptured", which was true of the file and false of the model.

Reads:  cplot_EQquasi.txt, or cplot_ft<N>_EQquasi.txt for --fault > 0
Writes: rupture_time.png per cycle
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import plotutils as pu
from seasio import read_array

pu.apply_style()
import matplotlib.pyplot as plt
import numpy as np

# .txt from a run, .csv if a gold directory ever carries one.
PATTERNS = ["cplot_EQquasi.*"]
UNRUPTURED = -1000.0
FNFT_COL = 15          # zero-based


def cplot_path(rdir, ift):
    """The cplot file for fault index `ift` (0-based)."""
    tag = "" if ift == 0 else f"ft{ift + 1}_"
    return os.path.join(rdir, f"cplot_{tag}EQquasi.txt")


def process_dir(rdir, outdir, interval, ift=0):
    path = cplot_path(rdir, ift)
    if not os.path.exists(path):
        print(f"  no {os.path.basename(path)} -- this run predates the "
              f"per-fault cplot (v1.13.0), or has fewer than {ift + 1} faults")
        return
    a = np.atleast_2d(read_array(path))
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

    # Column 2 of cplot_EQquasi.txt is -zcoor, i.e. depth below the free
    # surface, so larger means deeper. Invert the axis: drawn the natural way
    # up, the free surface is at the top and a rupture nucleating at 16 km and
    # growing towards 4 km reads as growing upwards, which is what it does.
    # cplot column 2 is -zcoor; negate to get the real coordinate.
    xk, zk = xs / 1e3, -zs / 1e3

    # Fixed contour interval, so the front's speed is directly readable: closely
    # spaced lines mean a slow front, widely spaced a fast one. A count-based
    # interval would change meaning between runs and make them incomparable.
    step = interval
    # Contour SECONDS SINCE THIS EVENT NUCLEATED, not absolute simulation
    # time. fnft is stamped on the cycle's own clock, so in a multi-cycle run
    # it carries the interseismic period too: cycle 2 of BP1002 nucleates at
    # 1.5e10 s, and 5-second contours drawn on that baseline all label as
    # "1.50757e+10" and overprint each other into a smear. Shifting also makes
    # cycles comparable, which absolute times never were.
    onset = float(np.nanmin(grid))
    grid = grid - onset
    lo, hi = float(np.nanmin(grid)), float(np.nanmax(grid))
    # A fixed interval keeps runs comparable, but draws nothing when the event
    # is shorter than one interval: BP7 ruptures over 1.1 s against BP5's 117,
    # so at 5 s the only levels are 0 and 5 and neither lies inside the data.
    # Fall back to a round interval giving ~8 lines, and say so on the figure.
    auto = False
    if (hi - lo) / step < 3.0:
        raw = max(hi - lo, 1e-12) / 8.0
        mag = 10.0 ** np.floor(np.log10(raw))
        step = next(m * mag for m in (1, 2, 2.5, 5, 10) if raw <= m * mag)
        auto = True
    clev = np.arange(np.floor(lo / step) * step, hi + step, step)

    # Labelled contour lines carry the times, so a colourbar repeats what the
    # labels already say and costs a fifth of the width. Drop it; keep a light
    # fill only to show the front's shape.
    fig, ax = plt.subplots(figsize=(8, 6), constrained_layout=True)
    cs = ax.contour(xk, zk, grid, levels=clev, colors="k", linewidths=1.1)
    ax.clabel(cs, inline=True, fontsize=11, fmt="%g")
    ax.set_xlabel("along strike (km)")
    ax.tick_params()
    # Column 2 of cplot_EQquasi.txt is -zcoor, so negate it back and plot the
    # real coordinate. Then the axis is just z, increasing upward, and nothing
    # here needs to know whether zero is a free surface (BP5, z in [-60, 0] km)
    # or the middle of a whole space (BP7 and BP8, z in [-0.5, 0.5] km). An
    # earlier version inverted the axis and annotated "free surface"
    # unconditionally, which drew BP7 upside down under a label that did not
    # apply.
    #
    # Full fault extent, not cropped to the ruptured patch: how much of the
    # fault did *not* rupture is part of the result.
    ax.set_xlim(xk.min(), xk.max())
    ax.set_ylim(zk.min(), zk.max())
    ax.set_ylabel("z (km)")
    ax.set_aspect("equal")
    ax.grid(alpha=0.15, lw=0.5)
    # The node count and time range are diagnostics, already printed to stdout;
    # in the title they crowd out the one thing a reader needs, which is what
    # the figure shows and which run it came from.
    ax.set_title(f"Rupture time since nucleation, contours every {step:g} s"
                 f"{' (auto: event shorter than the 5 s default)' if auto else ''}\n"
                 f"{os.path.basename(os.path.abspath(rdir))}"
                 f"  (nucleated at {onset:.6g} s)  fault {ift}")
    name = "rupture_time.png" if ift == 0 else f"rupture_time.f{ift}.png"
    pu.save(fig, pu.out_path(rdir, name, outdir), dpi=150)


def main():
    ap = pu.make_parser(__doc__, "plotRuptureTime.py",
                        "rupture_time.png into each cycle's directory "
                        "(empty run: a message, no figure).", PATTERNS)
    ap.add_argument("--fault", type=int, default=0,
                    help="fault index for multi-fault runs (default 0)")
    ap.add_argument("--interval", type=float, default=5.0, metavar="SECONDS",
                    help="contour interval in seconds (default 5). Fixed "
                         "rather than a level count, so line spacing means "
                         "the same thing across runs.")
    args = ap.parse_args()
    for label, rdir in pu.resolve_targets(args.dirs, PATTERNS,
                                          "plotRuptureTime.py"):
        print(f"processing {label or rdir} ({rdir})")
        process_dir(rdir, args.outdir, args.interval, args.fault)
    return 0


if __name__ == "__main__":
    sys.exit(main())
