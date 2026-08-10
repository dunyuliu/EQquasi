#! /usr/bin/env python3
"""Overlay the three readings of the BP8 initial condition.

BP8 prescribes tau_init (Table 1), V(0) = V_init (eq. 28) and
theta_0 = D_RS/V_init (eq. 30). Eq. (13) ties the three together, so only two
can hold. Each reading drops a different one, and the choice changes slip by
1.8x. This plots them side by side so the difference is visible rather than
tabulated.

Usage:
    plotICReadings.py [-o work/ic_readings.png]
"""

import argparse
import glob
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from seasio import read_rows, read_array

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# (directory, label, colour). Directories are gitignored run outputs.
READINGS = [
    ("work/c3.xi0.2",    "A literal: Table 1 tau + eq.(30) theta\n(violates eq.28: V0 = 65x V_init)", "tab:red"),
    ("work/bp8.icfix50", "C derived: eq.(28) V + Table 1 tau\n(theta_0 = 4.02e11 s)",                 "tab:orange"),
    ("work/bp8.equil50", "B equilibrated: eq.(28) V + eq.(30) theta\n(tau_0 = 12.9277 MPa)",          "tab:blue"),
]

REFERENCE_SLIP_MM = 21.0   # taehoKim_ref, read from a plot -- indicative only


def read(path):
    """Numeric rows of a benchmark output file. See scripts/seasio.py."""
    return read_array(path)
def load(d):
    g = os.path.join(ROOT, d, "global.dat")
    hits = glob.glob(os.path.join(ROOT, d, "fltst_strk+000dp+000.*"))
    if not os.path.exists(g) or not hits:
        return None
    return read(g), read(hits[0])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-o", "--out", default=os.path.join(ROOT, "work", "ic_readings.png"))
    args = ap.parse_args()

    cases = []
    for d, lab, c in READINGS:
        r = load(d)
        if r is None:
            print(f"skipping {d}: no output yet")
            continue
        cases.append((d, lab, c, r[0], r[1]))
    if not cases:
        print("no runs available")
        return 1

    fig, ax = plt.subplots(2, 2, figsize=(13, 8.5))
    for d, lab, c, g, s in cases:
        td = s[:, 0] / 86400.0
        ax[0, 0].plot(td, s[:, 1] * 1000, color=c, lw=1.6, label=lab)
        ax[0, 1].plot(g[:, 0] / 86400.0, g[:, 1], color=c, lw=1.6, label=lab)
        ax[1, 0].plot(td, s[:, 5], color=c, lw=1.6, label=lab)
        ax[1, 1].plot(td, s[:, 10], color=c, lw=1.6, label=lab)

    ax[0, 0].axhline(REFERENCE_SLIP_MM, color="k", ls=":", lw=1.2,
                     label=f"taehoKim_ref ~{REFERENCE_SLIP_MM:.0f} mm (from a plot)")
    # t_off = 100 hours.
    for a in ax.flat:
        a.axvline(100 / 24.0, color="gray", ls="--", lw=0.9, alpha=0.7)
        a.grid(alpha=0.3)
        a.set_xlabel("time (days)")

    ax[0, 0].set_ylabel("slip_2 at centre (mm)")
    ax[0, 0].set_title("Slip at the centre station", fontsize=10)
    ax[0, 1].set_ylabel("log10 Vmax (m/s)")
    ax[0, 1].set_title("Maximum slip rate over the frictional domain", fontsize=10)
    ax[1, 0].set_ylabel("shear_stress_2 (MPa)")
    ax[1, 0].set_title("Shear stress at the centre station", fontsize=10)
    ax[1, 1].set_ylabel("log10 state (s)")
    ax[1, 1].set_title("State variable at the centre station", fontsize=10)

    ax[0, 0].legend(fontsize=7, loc="lower right")
    fig.suptitle(
        "SEAS BP8-QD-GS: the three readings of an over-determined initial condition\n"
        "EQquasi, dx = 50 m, xi = 0.2. Dashed line: t_off = 4.17 d.",
        fontsize=11)
    fig.tight_layout()
    fig.savefig(args.out, dpi=115)
    print(f"wrote {args.out}")

    print()
    for d, lab, c, g, s in cases:
        print(f"{d:22s} slip={s[-1,1]*1000:7.2f} mm  peakVmax={g[:,1].max():.3f} "
              f"@ {g[g[:,1].argmax(),0]/86400:5.2f} d  V1={10**s[1,3]:.3e} m/s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
