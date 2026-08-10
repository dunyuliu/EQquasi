#! /usr/bin/env python3
"""Plot a run against its frozen gold reference, and show the residual.

`tests/e2e` tells you *whether* a run still matches its oracle. When it does not,
this tells you *where*. That gap is worth closing: a numeric assertion failure
names one station and one quantity, which is rarely enough to see whether a
change shifted a curve uniformly, moved a peak in time, or broke one region of
the fault.

Two layouts, matching the two kinds of oracle under reference/:

  BP5, BP5-dip90, BP7   an exact field snapshot at step 101. Plots run, gold and
                        (run - gold) as three maps over the fault plane, per
                        variable. Any non-zero residual is a regression -- these
                        references are compared at max-abs-diff 0.
  BP8                   station time series. Overlays run and gold per quantity
                        with the residual beneath, and prints the numbers the
                        e2e test asserts on.

Usage:
    plotAgainstGold.py bp5  work/mycase            [-o out_prefix]
    plotAgainstGold.py bp8  work/bp8.sub147        [-o out_prefix]
    plotAgainstGold.py bp7  work/mycase -v slip_rate

Exits non-zero if the run and the gold disagree beyond tolerance, so it can be
used as a check as well as a picture.
"""

import argparse
import glob
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from seasio import read_array


def read_gold_csv(path):
    """Gold references are comma-separated with a one-line header;
    run output is whitespace-separated with '#' comments. Different
    readers on purpose -- do not feed one to the other."""
    return np.genfromtxt(path, delimiter=',', skip_header=1)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

FIELD_BENCHMARKS = ("bp5", "bp5.dip90", "bp7")
SNAPSHOT = "fault.00101.nc"

# Station quantities worth eyeballing, as (column index, label).
BP8_PANELS = [(1, "slip_2 (m)"), (3, "slip_rate_2 (log10 m/s)"),
              (5, "shear_stress_2 (MPa)"), (7, "pore_pressure (MPa)"),
              (10, "state (log10 s)")]


def gold_dir(bench):
    d = os.path.join(ROOT, "reference", bench, "gold")
    if not os.path.isdir(d):
        raise SystemExit(f"no gold for '{bench}' at {d}")
    return d


def load_field(path):
    """Fault-plane snapshot as {variable: 2-D array}, from netCDF."""
    try:
        import netCDF4 as nc
    except ImportError:
        raise SystemExit("netCDF4 is needed to compare field snapshots")
    d = nc.Dataset(path)
    return {v: d.variables[v][:] for v in d.variables if not v.startswith("nid_")}


def compare_field(bench, run_dir, out, only=None):
    gold = load_field(os.path.join(gold_dir(bench), SNAPSHOT))
    hits = glob.glob(os.path.join(run_dir, "**", SNAPSHOT), recursive=True)
    if not hits:
        raise SystemExit(f"no {SNAPSHOT} under {run_dir}")
    run = load_field(hits[0])

    names = [only] if only else [v for v in gold if v in run]
    bad = []
    fig, ax = plt.subplots(len(names), 3, figsize=(13, 3.1 * len(names)),
                           squeeze=False)
    for r, v in enumerate(names):
        g, u = np.asarray(gold[v], float), np.asarray(run[v], float)
        if g.shape != u.shape:
            bad.append(f"{v}: shape {u.shape} vs gold {g.shape}")
            continue
        d = u - g
        if np.max(np.abs(d)) > 0.0:
            bad.append(f"{v}: max|diff| = {np.max(np.abs(d)):.3e}")
        for c, (arr, ttl, cmap) in enumerate(
                [(g, "gold", "viridis"), (u, "run", "viridis"),
                 (d, "run - gold", "RdBu_r")]):
            lim = np.max(np.abs(d)) or 1.0
            kw = dict(cmap=cmap, vmin=-lim, vmax=lim) if c == 2 else dict(cmap=cmap)
            m = ax[r, c].pcolormesh(arr, shading="nearest", **kw)
            plt.colorbar(m, ax=ax[r, c])
            ax[r, c].set_title(f"{v} — {ttl}", fontsize=9)
    fig.suptitle(f"{bench}: {run_dir} against reference/{bench}/gold", fontsize=11)
    fig.tight_layout()
    fig.savefig(f"{out}_field.png", dpi=110)
    print(f"wrote {out}_field.png")
    return bad


def compare_bp8(run_dir, out):
    g = gold_dir("bp8")
    summary = json.load(open(os.path.join(g, "summary.json")))
    stations = [k for k in summary if k.startswith(("+", "-"))]

    bad = []
    fig, ax = plt.subplots(len(BP8_PANELS), len(stations),
                           figsize=(6 * len(stations), 2.6 * len(BP8_PANELS)),
                           squeeze=False)
    for c, st in enumerate(stations):
        gold = read_gold_csv(os.path.join(g, f"fltst_strk{st}.csv"))
        hits = glob.glob(os.path.join(run_dir, f"fltst_strk{st}.*"))
        if not hits:
            bad.append(f"{st}: missing from the run")
            continue
        run = read_array(hits[0])
        for r, (col, label) in enumerate(BP8_PANELS):
            a = ax[r, c]
            a.plot(gold[:, 0] / 86400, gold[:, col], lw=2.4, alpha=.45, label="gold")
            a.plot(run[:, 0] / 86400, run[:, col], lw=1.2, label="run")
            a.set_ylabel(label, fontsize=8)
            a.grid(alpha=.3)
            if r == 0:
                a.set_title(f"station {st}", fontsize=10)
                a.legend(fontsize=7)
            if r == len(BP8_PANELS) - 1:
                a.set_xlabel("time (days)")
        # the two numbers the e2e test asserts on
        exp = summary[st]["slip_mm"]
        got = run[-1, 1] * 1000
        if abs(got - exp) / exp > 1e-4:
            bad.append(f"{st}: slip {got:.4f} mm vs gold {exp:.4f} mm")
    fig.suptitle(f"BP8: {run_dir} against reference/bp8/gold", fontsize=11)
    fig.tight_layout()
    fig.savefig(f"{out}_bp8.png", dpi=110)
    print(f"wrote {out}_bp8.png")
    return bad


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("benchmark", choices=list(FIELD_BENCHMARKS) + ["bp8"])
    ap.add_argument("run_dir")
    ap.add_argument("-o", "--out", default=None)
    ap.add_argument("-v", "--variable", default=None,
                    help="field benchmarks only: plot just this variable")
    args = ap.parse_args()
    out = args.out or os.path.join(args.run_dir, "vs_gold")

    if args.benchmark == "bp8":
        bad = compare_bp8(args.run_dir, out)
    else:
        bad = compare_field(args.benchmark, args.run_dir, out, args.variable)

    if bad:
        print("\nDISAGREES with gold:")
        for b in bad:
            print("  " + b)
        return 1
    print("\nmatches gold")
    return 0


if __name__ == "__main__":
    sys.exit(main())
