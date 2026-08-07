#! /usr/bin/env python3
"""Resample section 4.3 profiles onto the benchmark's 10 m output grid.

Section 4.3 requires profile data "with a spacing of 10 m (exactly)" from -400 m
to 400 m -- 81 nodes -- so that every participant's profiles land on a common
grid the platform can contour. That is a requirement on the *output sampling*,
not on the simulation resolution: section 6 explicitly invites results "from two
different spatial resolutions".

A run at dx = 50 m has 17 profile nodes and cannot satisfy section 4.3 directly.
This script linearly interpolates each time row onto the 81-node grid and writes
a submission-ready copy, leaving stations, global.dat and runInfo.json untouched.
Runs already on the 10 m grid are copied through unchanged.

The interpolation is recorded in each output header. It must be: the values are
resampled, not computed at those nodes, and a reader comparing our dx = 50 m
entry against a native dx = 10 m one is entitled to know that.

Usage:
    resampleBP8Profiles.py <result_dir> <out_dir>
"""

import argparse
import os
import shutil
import sys

import numpy as np

PROFILES = [
    "slip_2_strike.dat", "slip_2_depth.dat",
    "slip_3_strike.dat", "slip_3_depth.dat",
    "shear_stress_2_strike.dat", "shear_stress_2_depth.dat",
    "shear_stress_3_strike.dat", "shear_stress_3_depth.dat",
    "pore_pressure_strike.dat", "pore_pressure_depth.dat",
]

TARGET = np.arange(-400.0, 400.0 + 1e-9, 10.0)   # 81 nodes, exactly


def split_file(path):
    """Return (header_lines, fieldname_lines, numeric_rows)."""
    header, names, rows = [], [], []
    for line in open(path):
        s = line.rstrip("\n")
        t = s.strip()
        if not t:
            continue
        if t.startswith("#"):
            header.append(s)
            continue
        try:
            float(t.split()[0])
        except ValueError:
            names.append(s)
            continue
        rows.append([float(v) for v in t.split()])
    return header, names, rows


def resample(path, out_path):
    header, names, rows = split_file(path)
    d = np.array(rows)
    coords = d[0, 2:]
    if len(coords) == len(TARGET) and np.allclose(coords, TARGET):
        shutil.copy(path, out_path)
        return False

    out = np.zeros((d.shape[0], 2 + len(TARGET)))
    out[0, 0], out[0, 1] = 0.0, 0.0
    out[0, 2:] = TARGET
    for i in range(1, d.shape[0]):
        out[i, 0] = d[i, 0]          # time
        out[i, 1] = d[i, 1]          # max slip rate, taken over the whole fault
        out[i, 2:] = np.interp(TARGET, coords, d[i, 2:])

    with open(out_path, "w") as f:
        for h in header:
            # Correct the column-count line if it names the old width.
            if "Columns #3-" in h:
                f.write(f"# Columns #3-{2+len(TARGET):4d} = "
                        f"{h.split('=',1)[1].strip()}\n")
            else:
                f.write(h + "\n")
        f.write(f"# NOTE: profiles linearly interpolated from the "
                f"{len(coords)}-node simulation grid ({coords[1]-coords[0]:.0f} m "
                f"spacing) onto the {len(TARGET)}-node 10 m grid required by "
                f"section 4.3. Values are resampled, not computed at 10 m.\n")
        for n in names:
            f.write(n + "\n")
        f.write("# Here are the data\n")
        for i in range(out.shape[0]):
            f.write("".join(f"{v:22.14E}" for v in out[i]) + "\n")
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("result_dir")
    ap.add_argument("out_dir")
    args = ap.parse_args()

    if not os.path.isdir(args.result_dir):
        print(f"not a directory: {args.result_dir}")
        return 2
    os.makedirs(args.out_dir, exist_ok=True)

    # Everything that is not a profile passes through untouched.
    for name in sorted(os.listdir(args.result_dir)):
        src = os.path.join(args.result_dir, name)
        if not os.path.isfile(src):
            continue
        if name in PROFILES:
            continue
        if name.startswith("fltst_") or name in ("global.dat", "runInfo.json"):
            shutil.copy(src, os.path.join(args.out_dir, name))

    n = 0
    for name in PROFILES:
        src = os.path.join(args.result_dir, name)
        if not os.path.exists(src):
            print(f"missing {name}")
            return 1
        if resample(src, os.path.join(args.out_dir, name)):
            n += 1
    print(f"wrote {args.out_dir}: {n} profile(s) resampled onto the 81-node 10 m grid, "
          f"{len(PROFILES)-n} already on it")
    return 0


if __name__ == "__main__":
    sys.exit(main())
