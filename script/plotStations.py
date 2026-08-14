#! /usr/bin/env python3
"""Plot on-fault and off-fault station time series.

Every cycle writes two families of station file and until now nothing plotted
either of them:

BP5 and BP7 write these with a .txt extension, BP8 with .dat, so match on the
stem and accept either -- a single-extension glob silently finds nothing on the
other benchmark, which is indistinguishable from having no stations at all.

  fltst_strk<s>dp<d>.*        on-fault
  srfst_strk<s>st<b>dp<d>.*   off-fault, 7 columns

Columns, from `output_onfault_stations` / `output_offfault_stations`
(src/library_output.f90). The on-fault set differs between BP8 and the rest:
BP8 writes 11 columns with Darcy velocity and state, everything else writes 9.
The layout is detected from the column count, not from a flag.

  on fault (9)   t, slip_strike, -slip_dip, |V| (m/s), V_dip (m/s),
                 tau_strike (MPa), -tau_dip (MPa), sigma_eff (MPa),
                 log10 state.  Slip rates are LINEAR here.
  on fault (11)  t, slip_2, slip_3, log10 V_2, log10 V_3, tau_2, tau_3,
                 pore pressure, Darcy_2, Darcy_3, log10 state  [BP8]
  off fault (7)  t, h-disp, h-vel, v-disp, v-vel, n-disp, n-vel

Note slip_3 and tau_3 are written negated by the solver so that positive means
downwards; that sign convention is the file's, and is preserved here.

Usage:
    cd <case> && plotStations.py               all cycles
    cd <case> && plotStations.py cycle1        just that one
    plotStations.py <case>/cycle0             both families, one figure each
    plotStations.py <case>/cycle0 --off        off-fault only
"""

import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import plotutils as pu
from seasio import read_array

# One pattern, not three: resolve_results_dirs requires *all* patterns to
# match, and no run has both .txt and .dat station files -- BP5/BP7 write .txt,
# BP8 writes .dat. Both families end in "st_strk" (fltst_, srfst_).
PATTERNS = ("*st_strk*",)

# (column index, label) per layout. Time is column 0 throughout.
# 9-column layout (BP5, BP7). Read off the write statement in
# output_onfault_stations, NOT assumed from the BP8 one: the slip rates here
# are LINEAR m/s, column 7 is effective normal stress rather than pore
# pressure, and only the state variable is logged. Getting this wrong produces
# a figure that looks plausible and is mislabelled throughout.
#   t, slip_strike, -slip_dip, |V| (m/s), V_dip (m/s),
#   tau_strike (MPa), -tau_dip (MPa), sigma_eff (MPa), log10 state
ON_FAULT_9 = [(1, "slip, strike (m)"), (3, "slip rate (m/s)"),
              (5, r"$\tau$, strike (MPa)"), (7, r"$\bar\sigma$ (MPa)"),
              (8, r"log$_{10}\theta$ (s)")]
# 11-column layout (BP8, section 4.1). Slip rates ARE log10 here.
ON_FAULT_11 = [(1, "slip$_2$ (m)"), (3, r"log$_{10}$V$_2$ (m/s)"),
               (5, r"$\tau_2$ (MPa)"), (7, "pore pressure (MPa)"),
               (10, r"log$_{10}\theta$ (s)")]
OFF_FAULT = [(1, "horizontal disp (m)"), (2, "horizontal vel (m/s)"),
             (3, "vertical disp (m)"), (4, "vertical vel (m/s)"),
             (5, "normal disp (m)"), (6, "normal vel (m/s)")]


def station_files(rdir, off):
    """Every station file, including faults 2+.

    The pattern allows an optional ft<N>_ between the family prefix and
    'strk'. Faults 2 and beyond write fltst_ft2_strk..., and matching the
    literal prefix 'fltst_strk' skipped them entirely -- so on BP1002, which
    asks for three stations on each of two faults, this drew three and said
    nothing about the other three. A tool that silently halves a multi-fault
    model is worse than one that fails.
    """
    import re as _re
    family = "srfst" if off else "fltst"
    pat = _re.compile(rf"^{family}_(?:ft\d+_)?strk.*\.(?:txt|dat)$")
    hits = [f for f in sorted(os.listdir(rdir)) if pat.match(f)]
    return [os.path.join(rdir, f) for f in hits]


def panels_for(ncol, off):
    if off:
        return OFF_FAULT
    if ncol >= 11:
        return ON_FAULT_11
    return ON_FAULT_9


def station_label(path):
    """'strk000dp000' out of the filename -- the coordinates, without noise.

    An ft<N>_ tag is kept and moved to the front as 'fN ', because station
    coordinates are fault-local: fault 1 and fault 2 both have a strk000dp010,
    and two identically-labelled curves on one axes are worse than no label.
    """
    b = os.path.basename(path)
    for ext in (".txt", ".dat"):
        b = b.replace(ext, "")
    b = b.replace("fltst_", "").replace("srfst_", "")
    if b.startswith("ft"):
        tag, _, rest = b.partition("_")
        return f"f{tag[2:]} {rest}"
    return b


def plot(rdir, outdir, off, tunit):
    files = station_files(rdir, off)
    kind = "off-fault" if off else "on-fault"
    if not files:
        print(f"  no {kind} station files in {rdir}")
        return

    data = [(station_label(f), read_array(f)) for f in files]
    ncol = min(d.shape[1] for _, d in data)
    panels = [p for p in panels_for(ncol, off) if p[0] < ncol]
    div = {"s": 1.0, "d": 86400.0, "y": 86400.0 * 365.25}[tunit]
    xlab = {"s": "time (s)", "d": "time (days)", "y": "time (years)"}[tunit]

    # Enough height per panel that the curves are readable rather than
    # squeezed; line widths and fonts come from plotutils' shared style.
    fig, ax = plt.subplots(len(panels), 1, figsize=(9, 2.9 * len(panels)),
                           sharex=True, squeeze=False, constrained_layout=True)
    # Nine stations need a colour cycle that stays distinguishable; the default
    # ten-colour cycle repeats hues that are hard to tell apart in a line plot.
    colours = plt.get_cmap("turbo")(np.linspace(0.05, 0.95, len(data)))
    for r, (col, lab) in enumerate(panels):
        a = ax[r, 0]
        for (name, d), c in zip(data, colours):
            a.plot(d[:, 0] / div, d[:, col], lw=2.0, color=c, label=name)
        a.set_ylabel(lab)
        a.grid(alpha=0.25, lw=0.8)
        a.tick_params(length=5)
    ax[-1, 0].set_xlabel(xlab)

    # One legend for the figure: repeating nine station names in every panel
    # costs more space than the curves.
    ncols_leg = 3 if len(data) > 4 else 1
    # suptitle first, then the legend below it. Both at "outside upper center"
    # draw on top of each other -- which is what the first version did.
    fig.suptitle(f"{kind.capitalize()} stations — "
                 f"{os.path.basename(os.path.abspath(rdir))}")
    fig.legend(*ax[0, 0].get_legend_handles_labels(), loc="outside lower center",
               ncol=ncols_leg, frameon=False)

    name = ("off_fault_stations.png" if off else "on_fault_stations.png")
    pu.save(fig, pu.out_path(rdir, name, outdir))
    print(f"  {len(data)} {kind} stations, {ncol} columns")


def main():
    ap = pu.make_parser(__doc__, "plotStations.py",
                        "on_fault_stations.png / off_fault_stations.png into "
                        "each cycle's directory.", PATTERNS)
    ap.add_argument("--on", action="store_true",
                    help="only the on-fault (fltst_*) stations")
    ap.add_argument("--off", action="store_true",
                    help="only the off-fault (srfst_*) stations")
    ap.add_argument("--time", choices=("s", "d", "y"), default="d",
                    help="time axis unit: seconds, days (default) or years")
    args = ap.parse_args()
    # Both families by default -- a run has on-fault and off-fault stations and
    # you almost always want to see both. --on / --off narrow it.
    kinds = [False] if args.on else [True] if args.off else [False, True]
    for label, rdir in pu.resolve_targets(args.dirs, PATTERNS,
                                          "plotStations.py"):
        print(f"processing {label} ({rdir})")
        for off in kinds:
            plot(rdir, args.outdir, off, args.time)
    return 0


if __name__ == "__main__":
    sys.exit(main())
