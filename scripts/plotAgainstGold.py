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
import re
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

FIELD_BENCHMARKS = ("bp5", "bp5.dip90", "bp7", "stepover")
SNAPSHOT = "fault.00101.nc"

# Station quantities worth eyeballing, as (column index, label).
BP8_PANELS = [(1, "slip_2 (m)"), (3, "slip_rate_2 (log10 m/s)"),
              (5, "shear_stress_2 (MPa)"), (7, "pore_pressure (MPa)"),
              (10, "state (log10 s)")]


def gold_dir(bench):
    d = os.path.join(ROOT, "reference", bench)
    if not os.path.isdir(d):
        raise SystemExit(f"no gold for '{bench}' at {d}")
    return d


def load_field(path):
    """Fault-plane snapshot as {(variable, fault index): 2-D array}, from netCDF.

    netcdf_write_on_fault (src/netcdf_io.f90) always carries a leading
    nid_fault dimension of size ntotft, so single- and multi-fault runs write
    the same shape. Key on (name, ift) rather than squeezing that dimension
    away: a multi-fault snapshot then compares fault by fault, and a swapped
    or unwritten fault shows up as a named mismatch instead of hiding inside
    an aggregate. Single-fault cases just get ift = 0 throughout.
    """
    try:
        import netCDF4 as nc
    except ImportError:
        raise SystemExit("netCDF4 is needed to compare field snapshots")
    d = nc.Dataset(path)
    out = {}
    for v in d.variables:
        if v.startswith("nid_"):
            continue
        arr = d.variables[v][:]
        if arr.ndim == 3:
            for ift in range(arr.shape[0]):
                out[(v, ift)] = arr[ift]
        else:
            out[(v, 0)] = arr
    return out


def label(key):
    v, ift = key
    return v if ift == 0 else f"{v} [fault {ift}]"


# Units worth showing in, keyed by the netCDF variable name. Stresses in Pa
# print as "1e7" offsets that hide the value; slip rates and slip are already
# in sensible units.
_UNITS = {"shear_strike": (1e-6, "MPa"), "shear_dip": (1e-6, "MPa"),
          "effective_normal": (1e-6, "MPa"), "state_normal": (1e-6, "MPa"),
          "slip_rate": (1.0, "m/s"), "state_variable": (1.0, "s"),
          "slips": (1.0, "m"), "slipd": (1.0, "m"), "slipn": (1.0, "m")}


# Stress is conventionally shown on a blue-white-red diverging map. Applied to
# the stress fields and to residuals; sequential fields (slip rate, slip, state)
# keep a perceptually uniform map, where a diverging one would imply a
# meaningless centre.
_STRESS = ("shear_strike", "shear_dip", "effective_normal", "state_normal")


def _cmap_for(name, vmin, vmax):
    """(cmap, vmin, vmax). Stress is centred so white is the reference level:
    zero when the field changes sign, otherwise the field's own midpoint."""
    if name in _STRESS:
        if vmin < 0.0 < vmax:
            lim = max(abs(vmin), abs(vmax))
            return "bwr", -lim, lim
        mid = 0.5 * (vmin + vmax)
        half = max(vmax - mid, mid - vmin)
        return "bwr", mid - half, mid + half
    return "magma", vmin, vmax


def _scale(name, arr):
    return _UNITS.get(name, (1.0, ""))


def _case_params(run_dir):
    """Fault extent in metres from the run's model.txt, or None.

    model.txt is a positional contract (src/read_input.f90): line 1 is
    xmin xmax, line 3 is zmin zmax. Read it rather than importing the compset,
    so this works on any run directory including a bare cycle folder.
    """
    for d in (run_dir, os.path.dirname(os.path.normpath(run_dir))):
        f = os.path.join(d, "model.txt")
        if os.path.exists(f):
            try:
                rows = [l.split() for l in open(f) if l.strip()]
                return {"xlo": float(rows[0][0]), "xhi": float(rows[0][1]),
                        "zlo": float(rows[2][0]), "zhi": float(rows[2][1])}
            except (IndexError, ValueError):
                return None
    return None


def raw_dims(path):
    """Dimension sizes exactly as stored, before any per-fault reshaping."""
    import netCDF4 as nc
    d = nc.Dataset(path)
    return {k: len(v) for k, v in d.dimensions.items()}


def find_snapshot(run_dir):
    """(path, label) for the snapshot to compare, and say which cycle it is.

    A case run through run.sh keeps each cycle in Q<i>/, so the snapshot is
    almost never in the case root. This used to be a recursive glob taking
    hits[0] -- with Q0 and Q1 both present that silently compared whichever
    cycle the filesystem happened to return first, and reported "matches gold"
    or not depending on directory order. Resolve explicitly and name the cycle
    in the output instead.
    """
    direct = os.path.join(run_dir, SNAPSHOT)
    if os.path.exists(direct):
        return direct, os.path.basename(os.path.normpath(run_dir))

    cycles = sorted((d for d in glob.glob(os.path.join(run_dir, "cycle[0-9]*")) +
                     glob.glob(os.path.join(run_dir, "Q[0-9]*"))
                     if os.path.exists(os.path.join(d, SNAPSHOT))),
                    key=lambda d: int(re.sub(r"\D", "", os.path.basename(d)) or 0))
    if not cycles:
        raise SystemExit(
            f"no {SNAPSHOT} in {run_dir} or in any cycle*/Q* directory beneath it")
    if len(cycles) > 1:
        print(f"note: {len(cycles)} cycles have {SNAPSHOT} "
              f"({', '.join(os.path.basename(c) for c in cycles)}); "
              f"comparing {os.path.basename(cycles[0])}. "
              f"Pass the cycle directory to choose another.")
    return os.path.join(cycles[0], SNAPSHOT), os.path.basename(cycles[0])


def compare_field(bench, run_dir, out, only=None):
    gold_path = os.path.join(gold_dir(bench), SNAPSHOT)
    snap, cycle = find_snapshot(run_dir)
    hits = [snap]
    print(f"comparing {snap}")

    # Compare the stored shape BEFORE load_field reduces it. load_field keys on
    # (variable, fault index), which maps a 2-D (nz, nx) array and a 3-D
    # (1, nz, nx) array onto the same key -- so a gold file predating the
    # nid_fault dimension compared clean against a run that has it, and this
    # script reported "matches gold" while check.test.py's xarray identical()
    # correctly failed in CI. A comparator must not normalise away the
    # structural change it exists to detect.
    gd, rd = raw_dims(gold_path), raw_dims(hits[0])
    if gd != rd:
        raise SystemExit(
            f"{bench}: stored dimensions differ from gold, so the files are not "
            f"comparable.\n  gold {gold_path}: {gd}\n  run  {hits[0]}: {rd}\n"
            "Regenerate the gold from a current run.")

    gold = load_field(gold_path)
    run = load_field(hits[0])

    keys = [k for k in gold if k in run and (only is None or k[0] == only)]
    if only and not keys:
        raise SystemExit(f"no variable '{only}' in both run and gold")

    # A fault present in gold but missing from the run is a silent failure
    # mode worth naming: netcdf_write_on_fault sizes its output from ntotft,
    # so a run that meshed fewer faults writes fewer slabs.
    bad = [f"{label(k)}: missing from the run" for k in gold if k not in run]

    # Physical axes. The fault plane is tens of km across; drawn on array
    # indices with no aspect ratio it is unreadable, and a rectangular
    # nucleation patch reads as a diagonal smear.
    par = _case_params(run_dir)
    nz, nx = np.asarray(gold[keys[0]]).shape
    xs = (np.linspace(par["xlo"], par["xhi"], nx) / 1e3 if par else np.arange(nx))
    zs = (np.linspace(par["zlo"], par["zhi"], nz) / 1e3 if par else np.arange(nz))
    xlab, zlab = ("along strike (km)", "depth (km)") if par else ("strike index", "dip index")

    # Only draw a residual column when there is a residual. These oracles are
    # compared at max abs diff 0, so the common case is fifteen empty panels
    # taking a third of the figure to say "identical" -- state that in one line
    # instead. Gold and run share a colour scale, so they share one colourbar.
    resid = {k: float(np.max(np.abs(np.asarray(run[k], float)
                                    - np.asarray(gold[k], float))))
             for k in keys if np.asarray(gold[k]).shape == np.asarray(run[k]).shape}
    ncol = 2   # gold and run only; the residual is reported in the title and
               # on stdout, where a number is more use than a map of zeros.

    fig, ax = plt.subplots(len(keys), ncol,
                           figsize=(4.6 * ncol + 1.4, 2.6 * len(keys)),
                           squeeze=False, constrained_layout=True)
    for r, k in enumerate(keys):
        g, u = np.asarray(gold[k], float), np.asarray(run[k], float)
        if g.shape != u.shape:
            bad.append(f"{label(k)}: shape {u.shape} vs gold {g.shape}")
            for c in range(ncol):
                ax[r, c].set_visible(False)
            continue
        d, dmax = u - g, resid[k]
        if dmax > 0.0:
            bad.append(f"{label(k)}: max|diff| = {dmax:.3e}")
        sc, unit = _scale(k[0], g)
        vmin, vmax = float(min(g.min(), u.min())) * sc, float(max(g.max(), u.max())) * sc
        cmap, vmin, vmax = _cmap_for(k[0], vmin, vmax)
        for c, (arr, ttl) in enumerate([(g * sc, "gold"), (u * sc, "run")]):
            a_ = ax[r, c]
            m = a_.pcolormesh(xs, zs, arr, shading="nearest", cmap=cmap,
                              vmin=vmin, vmax=vmax)
            a_.set_aspect("equal")
            a_.set_xlabel(xlab, fontsize=8)
            a_.set_ylabel(zlab if c == 0 else "", fontsize=8)
            a_.tick_params(labelsize=7)
            a_.set_title(f"{label(k)} — {ttl}", fontsize=10)
        cb = fig.colorbar(m, ax=list(ax[r, :ncol]), fraction=0.025, pad=0.01)
        cb.set_label(unit or k[0], fontsize=8)
        cb.ax.tick_params(labelsize=7)
    worst = max(resid.values()) if resid else 0.0
    verdict = ("identical to gold — max|diff| = 0 on every variable"
               if worst == 0.0 else f"worst max|diff| = {worst:.3e}")
    fig.suptitle(f"{bench}  ·  {os.path.basename(os.path.normpath(run_dir))} "
                 f"{cycle}  against reference/{bench}/gold  ·  step 101, "
                 f"{nx} x {nz} fault nodes\n{verdict}", fontsize=12)
    fig.savefig(f"{out}_field.png", dpi=150)
    plt.close(fig)
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
    fig.suptitle(f"BP8: {run_dir} against reference/bp8", fontsize=11)
    fig.tight_layout()
    fig.savefig(f"{out}_bp8.png", dpi=110)
    print(f"wrote {out}_bp8.png")
    return bad


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="examples:\n"
               "  cd <case> && plotAgainstGold.py bp5      compare this case\n"
               "  plotAgainstGold.py bp5 work/mycase       from anywhere\n"
               "The field snapshot (fault.00101.nc) is searched recursively, "
               "so Q* cycle folders are found automatically.")
    ap.add_argument("benchmark", choices=list(FIELD_BENCHMARKS) + ["bp8"],
                    help="which frozen gold under reference/ to compare "
                         "against")
    ap.add_argument("run_dir", nargs="?", default=".",
                    help="run/case directory (default: current directory)")
    ap.add_argument("-o", "--out", default=None,
                    help="output figure prefix (default: <run_dir>/vs_gold)")
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
