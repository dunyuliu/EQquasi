#! /usr/bin/env python3
"""Plot one on-fault variable along a profile of the fault plane, per cycle.

Extracts a single profile from a fault snapshot (default fault.r.nc, the state
at the end of the cycle) and overlays every selected cycle on one axes:

  vertical         the column at mid-strike
  horizontal       the row at 10 km depth (measured down-dip; --dip corrects
                   the row spacing for dipping faults)
  sum_vertical /   the mean of all columns/rows whose mean value exceeds 0.5
  sum_horizontal   -- for 'slips', the average profile of the slipped region

For --var slips the field at nucleation is subtracted, so the profile shows
coseismic slip of this cycle's event: the snapshot closest to the nucleation
time (tdyna.txt) is used as reference. Without a tdyna.txt (no event), the
earliest snapshot is used and a note is printed.

Multi-fault runs: --fault selects which fault's plane is profiled (default 0).

Reads:  fault snapshot (--file), tdyna.txt, global.dat
Writes: profile.<var>.<type>.png (one figure overlaying all cycles)
"""

import glob
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import plotutils as pu

pu.apply_style()
import matplotlib.pyplot as plt
import numpy as np
import netCDF4 as nc

PATTERNS = ["fault.[0-9]*.nc", "global.dat"]


def slab_at_nucleation(rdir, var, ifault):
    """The variable's slab at the snapshot closest to nucleation."""
    snaps = [f for f in glob.glob(os.path.join(rdir, "fault.*.nc"))
             if pu.step_of(f) is not None]
    if not snaps:
        pu.die(f"plotProfile.py: no fault.NNNNN.nc snapshots in {rdir}")
    td = os.path.join(rdir, "tdyna.txt")
    ref = None
    if os.path.exists(td) and os.path.getsize(td) > 0:
        tdyna = np.atleast_1d(np.loadtxt(td))
        g = np.loadtxt(os.path.join(rdir, "global.dat"), ndmin=2)
        hit = np.where(g[:, 0] == tdyna[0])[0]
        if hit.size:
            ref = min(snaps, key=lambda f: abs(pu.step_of(f) - hit[0]))
    if ref is None:
        ref = min(snaps, key=pu.step_of)
        print(f"  no usable tdyna.txt in {rdir}; referencing slips to "
              f"{os.path.basename(ref)}")
    else:
        print(f"  slips referenced to nucleation-time snapshot "
              f"{os.path.basename(ref)}")
    with nc.Dataset(ref) as d:
        return pu.fault_slab(d.variables[var][:], ifault, var, ref)


def extract_profile(rdir, var, fileName, profile, ifault, dip, cellSize):
    path = os.path.join(rdir, fileName)
    if not os.path.exists(path):
        pu.die(f"plotProfile.py: no {fileName} in {rdir}")
    with nc.Dataset(path) as d:
        if var not in d.variables:
            pu.die(f"plotProfile.py: no variable '{var}' in {path}; "
                   f"available: {', '.join(d.variables)}")
        v = pu.fault_slab(d.variables[var][:], ifault, var, path)
    v = np.asarray(v, float)
    if var == "slips":
        v = v - np.asarray(slab_at_nucleation(rdir, var, ifault), float)

    if profile == "horizontal":
        step = cellSize / np.sin(np.radians(dip))
        row = v.shape[0] - round(10e3 / step)        # 10 km depth
        if not 0 <= row < v.shape[0]:
            pu.die(f"plotProfile.py: 10 km-depth row {row} outside the "
                   f"{v.shape[0]}-row grid; check --cellsize/--dip")
        return v[row, :]
    if profile == "vertical":
        return v[:, v.shape[1] // 2]
    axis, n = (0, v.shape[0]) if profile == "sum_horizontal" else (1, v.shape[1])
    lines = [np.take(v, i, axis=axis) for i in range(n)]
    lines = [l for l in lines if np.mean(l) > 0.5]
    return np.mean(lines, axis=0) if lines else \
        np.zeros(v.shape[1 - axis])


def main():
    ap = pu.make_parser(__doc__, "plotProfile.py",
                        "profile.<var>.<type>.png overlaying all selected "
                        "cycles (into the case dir, or cwd for a flat run).",
                        PATTERNS)
    ap.add_argument("--var", default="slips",
                    choices=["shear_strike", "shear_dip", "effective_normal",
                             "slips", "slipd", "slip_rate", "state_variable"],
                    help="on-fault variable to profile (default: slips, "
                         "shown as coseismic slip relative to nucleation)")
    ap.add_argument("--file", default="fault.r.nc",
                    help="snapshot to read (default fault.r.nc, end of cycle)")
    ap.add_argument("--profile_type", default="horizontal",
                    choices=["horizontal", "vertical",
                             "sum_horizontal", "sum_vertical"],
                    help="which profile to extract (default: horizontal, "
                         "at 10 km depth)")
    ap.add_argument("--fault", type=int, default=0,
                    help="fault index for multi-fault runs (default 0)")
    ap.add_argument("--dip", type=float, default=90.0,
                    help="fault dip in degrees, corrects the depth of the "
                         "horizontal profile (default 90)")
    ap.add_argument("--cellsize", type=float, default=2000.0,
                    help="on-fault cell size in metres (default 2000)")
    args = ap.parse_args()

    targets = pu.resolve_targets(args.dirs, PATTERNS, "plotProfile.py")
    scale, unit = 1.0, "m"
    if args.var in ("effective_normal", "shear_strike", "shear_dip"):
        scale, unit = 1e6, "MPa"
    elif args.var == "slip_rate":
        unit = "m/s"
    elif args.var == "state_variable":
        unit = "s"

    fig, ax = plt.subplots(figsize=(9, 5))
    for label, rdir in targets:
        print(f"processing {label or rdir} ({rdir})")
        p = extract_profile(rdir, args.var, args.file, args.profile_type,
                            args.fault, args.dip, args.cellsize)
        ax.plot(p / scale, lw=1.4,
                label=label or os.path.basename(os.path.abspath(rdir)))
    ax.set_xlabel(f"node # along {args.profile_type} profile")
    ax.set_ylabel(f"{args.var} ({unit})")
    ax.set_title(f"{args.var} — {args.profile_type} profile — "
                 f"{os.path.basename(args.file)}"
                 + (f" — fault {args.fault}" if args.fault else ""))
    ax.grid(alpha=0.3)
    ax.legend()
    fig.tight_layout()

    name = f"profile.{args.var}.{args.profile_type}.png"
    if args.outdir:
        os.makedirs(args.outdir, exist_ok=True)
        out = os.path.join(args.outdir, name)
    elif len(targets) == 1:
        out = os.path.join(targets[0][1], name)
    else:
        out = name
    pu.save(fig, out, dpi=150)
    return 0


if __name__ == "__main__":
    sys.exit(main())
