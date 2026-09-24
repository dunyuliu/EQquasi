#! /usr/bin/env python3
"""Event catalogue: moment magnitude vs time across cycles, one marker per event.

Each cycle is classified by which faults slipped seismically in it (any node
whose rupture time -- the first time its slip rate passed 1 mm/s -- was
stamped): one fault only, or JOINT when two or more faults did. A rupture that
jumps a step-over does so within one cycle, because a cycle only ends when slip
rate falls back after the event; a partner triggered days later lands in the
next cycle, and shows up in the lower panel as a short time since the previous
event.

    Mw = 2/3 (log10 M0 - 9.1),   M0 = mu * sum(area * coseismic slip)

over the nodes that ruptured, from cplot[_ftN]_ruptarea_trac_slip.txt, with
mu = rho * vs**2 from the case's user_defined_params.py. Output that predates
that file falls back to end-of-cycle slip on the ruptured nodes times dx*dz,
and says so. A fault whose per-fault rupture-time file is missing (output older
than v1.13) is reported as unknown, never as quiet.

Several directories may be given, in order, to join a run and its restarts
(e.g. a run and a fork continued from one of its cycles). A cycle number that
appears again in a later directory replaces the earlier one, which is how a
fork overrides the history it continued from.

Reads:  global.dat (cycle length; not needed for the last cycle), cplot[_ftN]_EQquasi.txt (rupture time),
        cplot[_ftN]_ruptarea_trac_slip.txt (area, slip)
Writes: mag_vs_time.png in the case root (above result/cycleN), or -o DIR
"""

import glob
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import plotutils as pu

pu.apply_style()
import matplotlib.pyplot as plt
import numpy as np

PATTERNS = ["cplot_EQquasi.*"]
UNRUPTURED = -1000.0   # plotRuptureTime.py's sentinel for "never ruptured"
FNFT_COL = 15          # zero-based rupture-time column of cplot_EQquasi
YEAR = 365.25 * 86400.0
TOOL = "plotMagnitudeTime.py"


def tag(ift):
    return "" if ift == 0 else f"ft{ift + 1}_"


def find(rdir, stem):
    hits = sorted(glob.glob(os.path.join(rdir, stem + ".*")))
    return hits[0] if hits else None


def cycle_number(rdir):
    m = re.search(r"(?:cycle|Q)(\d+)$", os.path.basename(os.path.abspath(rdir)))
    return int(m.group(1)) if m else -1


def cycle_events(rdir, par, mu, nft):
    """[(fault, onset s, M0, source note)] for one cycle, plus unknown faults."""
    ev, unknown = [], []
    for ift in range(nft):
        cp = find(rdir, f"cplot_{tag(ift)}EQquasi")
        if cp is None:
            unknown.append(ift)
            continue
        fnft = np.atleast_2d(np.loadtxt(cp))[:, FNFT_COL]
        hit = fnft > UNRUPTURED / 2.0
        if not hit.any():
            continue
        if mu is None:
            pu.die(f"{TOOL}: {rdir} has a seismic event but no "
                   "user_defined_params.py to take mu from; pass --mu.")
        ra = find(rdir, f"cplot_{tag(ift)}ruptarea_trac_slip")
        if ra is not None:
            a = np.atleast_2d(np.loadtxt(ra))
            m0, note = mu * np.sum(a[:, 0] * a[:, 1]), ""
        else:
            if par is None:
                pu.die(f"{TOOL}: {rdir} has no ruptarea file and no case "
                       "parameters for the node area; cannot compute a moment.")
            m0, note = mu * end_slip(rdir, ift, hit) * par.dx * par.dz, \
                "end-of-cycle slip (no ruptarea file)"
        ev.append((ift, fnft[hit].min(), m0, note))
    return ev, unknown


def end_slip(rdir, ift, hit):
    """Sum of end-of-cycle slip over the ruptured nodes of fault ift."""
    import xarray as xr
    rnc = os.path.join(rdir, "fault.r.nc")
    if not os.path.exists(rnc):
        pu.die(f"{TOOL}: {rdir} has neither a ruptarea file nor fault.r.nc; "
               "cannot compute a moment.")
    s = np.abs(xr.open_dataset(rnc)["slips"].values)
    s = s[ift] if s.ndim == 3 else s
    flat = s.reshape(-1)
    return float(np.sum(flat[:hit.size][hit])) if flat.size >= hit.size else 0.0


def main():
    ap = pu.make_parser(__doc__, TOOL,
                        "mag_vs_time.png in the case root (above result/cycleN).",
                        PATTERNS)
    ap.add_argument("--mu", type=float, default=None,
                    help="shear modulus, Pa (default: rho*vs**2 from the case)")
    ap.add_argument("--pair-yr", type=float, default=5.0,
                    help="lower-panel band for a triggered partner, yr "
                         "(default 5)")
    args = ap.parse_args()

    targets = pu.resolve_targets(args.dirs, PATTERNS, TOOL)
    by_cycle = {}
    for _, rdir in targets:                 # later directories override
        by_cycle[cycle_number(rdir)] = rdir
    cycles = sorted(by_cycle)
    # Case parameters give mu and the fault count. A frozen gold directory
    # may carry no user_defined_params.py; then the fault count comes from the
    # per-fault files, and mu is demanded only if there is a moment to compute.
    d0 = os.path.abspath(by_cycle[cycles[0]])
    here = [d0, os.path.dirname(d0), os.path.dirname(os.path.dirname(d0))]
    par = (pu.load_par(d0, TOOL)
           if any(os.path.isfile(os.path.join(d, "user_defined_params.py"))
                  for d in here) else None)
    nft = (int(getattr(par, "ntotft", 1)) if par is not None else
           1 + len(glob.glob(os.path.join(d0, "cplot_ft[0-9]*_EQquasi.*"))))
    mu = (args.mu if args.mu is not None else
          par.rou * par.vs ** 2 if par is not None else None)

    rows, notes, t0 = [], set(), 0.0
    for c in cycles:
        rdir = by_cycle[c]
        # Cycle length only offsets the cycles after it, so a lone cycle
        # without global.dat (a flat BP8 gold directory) needs none.
        g = os.path.join(rdir, "global.dat")
        if os.path.exists(g):
            dur = float(np.atleast_2d(np.loadtxt(g))[-1, 0])
        elif c == cycles[-1]:
            dur = 0.0
        else:
            pu.die(f"{TOOL}: {rdir} has no global.dat, so the cycles after "
                   "it cannot be placed in time.")
        ev, unknown = cycle_events(rdir, par, mu, nft)
        for u in unknown:
            notes.add(f"fault {u + 1}: no per-fault rupture-time file "
                      "(output predates v1.13), shown as unknown")
        for e in ev:
            if e[3]:
                notes.add(e[3])
        if ev:
            faults = sorted({e[0] for e in ev})
            m0 = sum(e[2] for e in ev)
            kind = "joint" if len(faults) > 1 else f"fault {faults[0] + 1}"
            rows.append((c, kind, (t0 + min(e[1] for e in ev)) / YEAR,
                         (2.0 / 3.0) * (np.log10(m0) - 9.1) if m0 > 0 else np.nan))
        t0 += dur

    kinds = [f"fault {i + 1}" for i in range(nft)] + ["joint"]
    colour = {k: f"C{i}" for i, k in enumerate(kinds[:-1])}
    colour["joint"] = "C3"
    print(f"{len(cycles)} cycles, {len(rows)} events: " +
          ", ".join(f"{k} {sum(r[1] == k for r in rows)}" for k in kinds))
    for n in sorted(notes):
        print("  note: " + n)

    fig, (ax, ax2) = plt.subplots(2, 1, figsize=(10, 5.4), sharex=True,
                                  gridspec_kw={"height_ratios": [1.2, 1]})
    for k in kinds:
        e = [r for r in rows if r[1] == k]
        t, m = [r[2] for r in e], [r[3] for r in e]
        if e:
            lo = np.nanmin([r[3] for r in rows]) - 0.1
            ax.vlines(t, lo, m, color=colour[k], lw=1.1)
        ax.scatter(t, m, color=colour[k], s=26, zorder=3,
                   label=f"{k} ({len(e)})")
    for i in range(1, len(rows)):
        dt = rows[i][2] - rows[i - 1][2]
        ax2.scatter(rows[i][2], dt, color=colour[rows[i][1]], s=26, zorder=3)
    if rows:
        ax2.set_yscale("log")
        ax2.axhspan(ax2.get_ylim()[0], args.pair_yr, color="0.9", zorder=0)
    else:
        ax.text(0.5, 0.5, "no seismic events (no node passed 1 mm/s)",
                transform=ax.transAxes, ha="center")
    ax.set_ylabel("Mw")
    ax.legend(fontsize=7, loc="lower right", ncol=len(kinds))
    ax2.set_ylabel("time since previous event (yr)")
    ax2.set_xlabel("time (yr)")
    for a in (ax, ax2):
        a.grid(alpha=0.3, which="both")
    title = f"{len(cycles)} cycles, {len(rows)} events"
    if notes:
        title += "  [" + "; ".join(sorted(notes)) + "]"
    ax.set_title(title, fontsize=8)

    # A multi-cycle summary belongs to the case, not to one cycle: write it
    # in the case root (above result/cycleN), or -o DIR.
    first = os.path.abspath(by_cycle[cycles[0]])
    home = os.path.dirname(first) if cycle_number(first) >= 0 else first
    if os.path.basename(home) == "result":
        home = os.path.dirname(home)
    out = (os.path.join(home, "mag_vs_time.png") if args.outdir is None
           else pu.out_path(home, "mag_vs_time.png", args.outdir))
    pu.save(fig, out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
