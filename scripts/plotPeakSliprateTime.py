#! /usr/bin/env python3
"""Peak slip rate against time, concatenated across earthquake cycles.

Reads global.dat from each cycle (Q0, Q1, ... in numerical order -- cycle
time restarts at zero each cycle, so times are accumulated to a single
continuous axis) and plots max slip rate on a log axis against time in years.
Naming several cases overlays them, one line per case, for comparing e.g.
different dips or resolutions.

Column 2 of global.dat is the linear max slip rate in m/s for every benchmark
except BP8, which stores log10(m/s) (src/library_output.f90); both are shown
on the same log axis by converting the BP8 form back to linear. The BP8 case
is detected from the file's own '# Column #2 ... log10' header, not guessed.

Reads:  global.dat per cycle
Writes: peak_slip_rate_vs_time.png (into the case directory; cwd when several
        cases are overlaid)
"""

import glob
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import plotutils as pu

pu.apply_style()
import matplotlib.pyplot as plt
import numpy as np

PATTERNS = ["global.dat"]
SECONDS_PER_YEAR = 365.25 * 24 * 3600


def load_case(tokens_or_dir):
    """(t_years, vmax) concatenated over the case's cycles."""
    parts = pu.resolve_targets(tokens_or_dir, PATTERNS, "plotPeakSliprateTime.py")
    times, rates, t0 = [], [], 0.0
    for label, rdir in parts:
        path = os.path.join(rdir, "global.dat")
        is_log10 = False
        with open(path) as f:
            for line in f:
                if line.startswith("#") and "Column #2" in line:
                    is_log10 = "log10" in line
        d = np.loadtxt(path, ndmin=2)
        v = 10.0 ** d[:, 1] if is_log10 else d[:, 1]
        times.append(d[:, 0] + t0)
        rates.append(v)
        t0 += d[-1, 0]
        print(f"  {label or rdir}: {len(d)} steps, "
              f"{d[-1, 0] / 86400:.2f} d simulated")
    return np.concatenate(times) / SECONDS_PER_YEAR, np.concatenate(rates)


def main():
    ap = pu.make_parser(__doc__, "plotPeakSliprateTime.py",
                        "peak_slip_rate_vs_time.png, one line per named case "
                        "with all its cycles concatenated.", PATTERNS)
    args = ap.parse_args()

    # Group tokens into cases: a case directory (or bare cycle names, which
    # all belong to the case in cwd) each becomes one line.
    tokens = args.dirs or [""]
    case_dirs = []
    bare = [t for t in tokens if t and not os.path.isdir(t)
            and os.path.isdir(f"Q{t.lstrip('Q')}")]
    if bare:                      # named cycles of the case in cwd: one line
        case_dirs.append((os.path.basename(os.getcwd()), bare))
        tokens = [t for t in tokens if t not in bare]
    for t in tokens:
        if t == "":
            case_dirs.append((os.path.basename(os.getcwd()), []))
        else:
            case_dirs.append((os.path.basename(os.path.abspath(t)), [t]))

    fig, ax = plt.subplots(figsize=(9, 5))
    for name, toks in case_dirs:
        t, v = load_case(toks)
        ax.plot(t, np.maximum(v, 1e-30), lw=1.2, label=name)
    ax.set_yscale("log")
    ax.set_xlabel("time (years)")
    ax.set_ylabel("peak slip rate (m/s)")
    ax.set_title("Peak slip rate vs time")
    ax.grid(alpha=0.3)
    ax.legend()
    fig.tight_layout()

    if args.outdir:
        out = os.path.join(args.outdir, "peak_slip_rate_vs_time.png")
        os.makedirs(args.outdir, exist_ok=True)
    elif len(case_dirs) == 1 and case_dirs[0][1] and \
            os.path.isdir(case_dirs[0][1][0]):
        out = os.path.join(case_dirs[0][1][0], "peak_slip_rate_vs_time.png")
    else:
        out = "peak_slip_rate_vs_time.png"
    pu.save(fig, out, dpi=150)
    return 0


if __name__ == "__main__":
    sys.exit(main())
