#! /usr/bin/env python3
"""Overlay the campaign's domain-size and xi cases on one set of axes.

Reads work/camp.dom*/ and work/camp.xi*/ and writes:
    work/campaign_domain.png   key variables vs time, one line per domain size
    work/campaign_xi.png       same for the time-step safety factor
"""

import glob
import os
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from seasio import read_rows, read_array

R = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def rd(path):
    """Numeric rows of a benchmark output file. See scripts/seasio.py."""
    return read_array(path)
def collect(pattern, keyfn):
    out = []
    for d in sorted(glob.glob(os.path.join(R, "work", pattern)), key=keyfn):
        g = os.path.join(d, "global.dat")
        # BP8 station files are .dat, not .txt -- the platform routes on the
        # extension and ignores .txt. This looked for .txt only, so it silently
        # matched nothing and reported "no cases" for every BP8 sweep, which is
        # indistinguishable from having no runs. Accept either.
        hits = glob.glob(os.path.join(d, "fltst_strk+000dp+000.*"))
        if os.path.exists(g) and hits:
            out.append((d, rd(g), rd(hits[0])))
    return out


def panel(cases, labelfn, title, outfile):
    if not cases:
        print(f"no cases for {outfile}")
        return
    fig, ax = plt.subplots(2, 2, figsize=(12, 8))
    for d, g, s in cases:
        lab = labelfn(d)
        ax[0, 0].plot(g[:, 0], g[:, 1], lw=1.5, label=lab)
        ax[0, 1].plot(g[:, 0], g[:, 2], lw=1.5, label=lab)
        ax[1, 0].plot(s[:, 0], s[:, 1] * 1000, lw=1.5, label=lab)
        ax[1, 1].plot(s[:, 0], s[:, 5], lw=1.5, label=lab)
    for a, t, yl in (
        (ax[0, 0], "Maximum slip rate", "log10 V (m/s)"),
        (ax[0, 1], "Moment density rate", "N/s"),
        (ax[1, 0], "Slip at centre station", "slip_2 (mm)"),
        (ax[1, 1], "Shear stress at centre station", "shear_stress_2 (MPa)"),
    ):
        a.set_title(t, fontsize=10)
        a.set_xlabel("time (s)")
        a.set_ylabel(yl, fontsize=9)
        a.grid(alpha=0.3)
        a.legend(fontsize=8)
    fig.suptitle(title, fontsize=11)
    fig.tight_layout()
    fig.savefig(outfile, dpi=115)
    print(f"wrote {outfile}")


DOM_RE = re.compile(r"dom(\d+)")
XI_RE = re.compile(r"xi([\d.]+)")


def dom_key(d):
    return int(DOM_RE.search(d).group(1))


def xi_key(d):
    return float(XI_RE.search(d).group(1))


def main():
    dom = collect("c*.dom*", dom_key)
    panel(dom, lambda d: f"half-width {dom_key(d)} m",
          "BP8-QD-GS, dx = 50 m: sensitivity to elastic domain size",
          os.path.join(R, "work", "campaign_domain.png"))

    xi = collect("c*.xi*", xi_key)
    panel(xi, lambda d: f"xi = {xi_key(d)}",
          "BP8-QD-GS, 500 m box: sensitivity to the time-step factor xi",
          os.path.join(R, "work", "campaign_xi.png"))


if __name__ == "__main__":
    main()
