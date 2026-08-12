#! /usr/bin/env python3
"""Shared helpers for the post-processing utilities in scripts/.

Three conventions every utility must follow live here, once, because each has
already been violated by copies drifting apart:

1. The on-fault netCDF contract. `netcdf_write_on_fault` (src/netcdf_io.f90)
   stores every on-fault variable 3-D as (nid_fault, nid_dip, nid_strike),
   with nid_fault of size ntotft even for single-fault runs. Files written
   before the mfault merge are 2-D (nid_dip, nid_strike). Every consumer goes
   through `fault_slab()` so that a multi-fault run is plotted fault by fault
   (never silently reduced to the first fault, or fed 3-D to a 2-D plot call,
   which xarray turns into a histogram), while legacy 2-D files keep working.

2. The directory argument. Utilities take a case or results directory
   (default '.') and read nothing implicitly from the current working
   directory. A case run through run.sh puts each cycle's output in Q<i>/
   subdirectories, with only the restart files (fault.r.nc, disp.r.nc) copied
   back to the case root -- so `resolve_results_dirs()` treats a directory as
   results if it holds the files the tool actually reads, otherwise discovers
   Q* cycles (numerically sorted; restarted cases may begin at Q3), and
   otherwise stops with the missing file named. `load_par()` resolves
   user_defined_params.py from the results directory or its parent (the case
   root, for a Q<i> cycle folder), never from sys.path luck.

3. Figure style. 150 dpi, labelled axes with units, labelled colourbars,
   a title naming case / fault / step / time, diverging colourmaps centred on
   the reference value for stress-like fields. `apply_style()` + `save()`
   keep that uniform without inventing a plotting framework.
"""

import glob
import os
import re
import sys

import numpy as np


# One place for figure style, applied on import so every utility that uses this
# module gets the same look. Sizes are chosen for figures printed at roughly
# 8 x 6 inches: readable in a talk or a paper column without rescaling, which
# the matplotlib defaults (10 pt on a 6.4 x 4.8 in figure) are not.
def apply_style():
    import matplotlib
    matplotlib.rcParams.update({
        "figure.dpi": 150,
        "savefig.dpi": 150,
        "font.size": 14,
        "axes.titlesize": 16,
        "axes.labelsize": 15,
        "xtick.labelsize": 13,
        "ytick.labelsize": 13,
        "legend.fontsize": 13,
        "figure.titlesize": 17,
        "axes.linewidth": 1.1,
        "lines.linewidth": 2.2,
        "lines.markersize": 8,
        "xtick.major.width": 1.1,
        "ytick.major.width": 1.1,
        "xtick.major.size": 5.5,
        "ytick.major.size": 5.5,
        "grid.alpha": 0.25,
        "savefig.bbox": "tight",
    })


apply_style()


def die(msg):
    """Fail with a message, not a traceback. These tools are run by people
    looking at results, not reading source."""
    raise SystemExit(msg)


# ---------------------------------------------------------------- fault slabs

def nfault(arr):
    """Number of faults stored in an on-fault variable (1 for legacy 2-D)."""
    a = np.asarray(arr)
    return a.shape[0] if a.ndim == 3 else 1


def fault_slab(arr, ifault=0, name="variable", path=""):
    """The (nid_dip, nid_strike) slab of one fault, from a 2-D or 3-D array."""
    a = np.asarray(arr)
    where = f"{path}: " if path else ""
    if a.ndim == 2:
        if ifault != 0:
            die(f"{where}{name} is 2-D (pre-nid_fault layout), so only fault 0 "
                f"exists; got fault index {ifault}")
        return a
    if a.ndim == 3:
        if not 0 <= ifault < a.shape[0]:
            die(f"{where}{name} holds {a.shape[0]} fault(s) along nid_fault; "
                f"fault index {ifault} is out of range")
        return a[ifault]
    die(f"{where}{name} has unexpected shape {a.shape}; expected "
        "(nid_dip, nid_strike) or (nid_fault, nid_dip, nid_strike)")


# ------------------------------------------------------- directory resolution

def resolve_results_dirs(path, patterns, tool=""):
    """[(cycle_label, dir)] of results directories reachable from `path`.

    `patterns` are the glob patterns of files the tool actually reads
    (e.g. ["global.dat", "fault.[0-9]*.nc"]). Requiring those -- not just any
    output -- matters because run.sh copies fault.r.nc/disp.r.nc back to the
    case root between cycles, which must not make the root look like results.
    """
    tag = f"{tool}: " if tool else ""
    path = path or "."
    if not os.path.isdir(path):
        die(f"{tag}not a directory: {path}")

    def ok(d):
        return all(glob.glob(os.path.join(d, p)) for p in patterns)

    if ok(path):
        return [("", path)]
    qs = []
    # run.sh writes cycle0, cycle1, ... Older cases used Q0, Q1 -- accept both
    # so existing output stays readable.
    for d in glob.glob(os.path.join(path, "cycle[0-9]*")) + \
             glob.glob(os.path.join(path, "Q[0-9]*")):
        m = re.fullmatch(r"(?:cycle|Q)(\d+)", os.path.basename(d))
        if m and os.path.isdir(d) and ok(d):
            qs.append((int(m.group(1)), d))
    if qs:
        return [(os.path.basename(d), d) for _, d in sorted(qs)]
    die(f"{tag}found none of the required files ({', '.join(patterns)}) in "
        f"{path} or in any {os.path.join(path, 'Q*')} cycle directory.\n"
        "Point me at a results directory, or at a case directory whose Q* "
        "cycle folders contain the run output.")


def resolve_targets(tokens, patterns, tool=""):
    """[(label, dir)] from zero or more user-supplied directory/cycle tokens.

    No tokens: discover everything from '.', per resolve_results_dirs.
    Each token may be a results directory, a case directory (expands to all
    its Q* cycles), a cycle name like 'Q1', or a bare number '1'. Naming
    cycles processes exactly those. A token that resolves to nothing is an
    error listing the cycles that do exist.
    """
    tag = f"{tool}: " if tool else ""
    if not tokens:
        return resolve_results_dirs(".", patterns, tool)
    out = []
    for t in tokens:
        cand = t
        if not os.path.isdir(cand):
            m = re.fullmatch(r"(?:cycle|Q)?(\d+)", t)
            hit = next((f"{pre}{int(m.group(1))}" for pre in ("cycle", "Q")
                        if m and os.path.isdir(f"{pre}{int(m.group(1))}")), None)
            if hit:
                cand = hit
            else:
                have = sorted(os.path.basename(d)
                              for d in glob.glob("cycle[0-9]*") + glob.glob("Q[0-9]*")
                              if os.path.isdir(d))
                die(f"{tag}no such directory or cycle: {t}"
                    + (f" (cycles here: {', '.join(have)})" if have else ""))
        out.extend(resolve_results_dirs(cand, patterns, tool))
    return out


def out_path(results_dir, filename, outdir=None):
    """Where a figure for `results_dir` goes: beside the input by default,
    into `outdir` (suffixed with the cycle/case name, so cycles stay
    distinguishable) when the user redirected output."""
    if outdir is None:
        return os.path.join(results_dir, filename)
    os.makedirs(outdir, exist_ok=True)
    base = os.path.basename(os.path.abspath(results_dir))
    root, ext = os.path.splitext(filename)
    return os.path.join(outdir, f"{root}_{base}{ext}")


def make_parser(doc, prog, writes, patterns):
    """The uniform CLI: DIR_OR_CYCLE..., -o, and an epilog that spells out the
    cycle behaviour -- `-h` is the only place a user standing in a case folder
    can learn that the bare command processes every cycle."""
    import argparse
    epilog = f"""\
cycle selection (a case run via run.sh keeps each cycle's output in Q<i>/):
  cd <case> && {prog}              all cycles (cycle0, cycle1, ...) in the case
  cd <case> && {prog} cycle1           just Q1 (bare '1' also works)
  cd <case> && {prog} cycle1 Q2        exactly those two
  {prog} <case>                    all cycles of a case, from anywhere
  {prog} <case>/cycle0            one cycle, from anywhere
A directory that itself holds the required files ({', '.join(patterns)})
is treated as a single results directory and processed directly.

output:
  {writes}
  With -o DIR everything goes there instead, suffixed with the cycle name."""
    p = argparse.ArgumentParser(
        prog=prog, description=doc,
        formatter_class=argparse.RawDescriptionHelpFormatter, epilog=epilog)
    p.add_argument("dirs", nargs="*", metavar="DIR_OR_CYCLE",
                   help="results/case directories or cycle names (Q1 or 1); "
                        "none = process every cycle found from the current "
                        "directory")
    p.add_argument("-o", "--outdir", default=None, metavar="DIR",
                   help="write output here instead of beside the input "
                        "(default: into each results directory)")
    return p


def load_par(run_dir, tool=""):
    """The `par` object from user_defined_params.py in run_dir or its parent.

    A Q<i> cycle folder holds output only; the parameters live one level up in
    the case root. Loads by explicit path (never via cwd or PYTHONPATH), with
    scripts/ temporarily on sys.path for `from defaultParameters import ...`.
    """
    import importlib.util
    tag = f"{tool}: " if tool else ""
    run_dir = os.path.abspath(run_dir or ".")
    for d in (run_dir, os.path.dirname(run_dir)):
        p = os.path.join(d, "user_defined_params.py")
        if os.path.isfile(p):
            break
    else:
        die(f"{tag}no user_defined_params.py in {run_dir} or its parent; "
            "this tool needs the case parameters (point it at the case or "
            "cycle directory).")
    added = [os.path.dirname(os.path.abspath(__file__)),
             os.path.dirname(p)]
    sys.path[:0] = added
    stale = sys.modules.pop("user_defined_params", None)
    dwb = sys.dont_write_bytecode
    sys.dont_write_bytecode = True
    try:
        spec = importlib.util.spec_from_file_location("user_defined_params", p)
        mod = importlib.util.module_from_spec(spec)
        try:
            spec.loader.exec_module(mod)
        except Exception as e:                       # noqa: BLE001
            die(f"{tag}failed to load {p}: {e}")
        return mod.par
    finally:
        sys.dont_write_bytecode = dwb
        sys.modules.pop("user_defined_params", None)
        if stale is not None:
            sys.modules["user_defined_params"] = stale
        for a in added:
            if a in sys.path:
                sys.path.remove(a)


def fault_coords(par, ndip, nstrike):
    """(x, z, xlabel, ylabel) for fault-plane axes, in km when the case's
    fault grid matches the snapshot, node indices otherwise (e.g. multi-fault
    cases, where par.fx spans the union of segments, not one fault)."""
    try:
        fx, fz = np.asarray(par.fx), np.asarray(par.fz)
        if fx.size == nstrike and fz.size == ndip:
            return fx / 1e3, fz / 1e3, "along strike (km)", "z (km)"
    except AttributeError:
        pass
    return (np.arange(nstrike), np.arange(ndip),
            "node # along strike", "node # along dip")


# ----------------------------------------------------------------- aesthetics

def apply_style():
    import matplotlib
    matplotlib.use("Agg")
    matplotlib.rcParams.update({
        "font.size": 10,
        "axes.titlesize": 11,
        "axes.labelsize": 10,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "legend.fontsize": 9,
        "figure.titlesize": 13,
        "axes.grid": False,
    })


def save(fig, path, dpi=150):
    """Save, close, and say where it went (with size: keep single panels
    under ~500 KB -- raise dpi only where a figure genuinely needs it)."""
    import matplotlib.pyplot as plt
    fig.savefig(path, dpi=dpi, pil_kwargs={"optimize": True})
    plt.close(fig)
    print(f"wrote {path} ({os.path.getsize(path) / 1024:.0f} KB)")


# -------------------------------------------------------------- odds and ends

def step_of(fname):
    """Time-step number encoded in fault.NNNNN.nc / disp.NNNNN.nc."""
    m = re.search(r"\.(\d+)\.nc$", os.path.basename(fname))
    return int(m.group(1)) if m else None


def _alphanum_key(s):
    return [int(c) if c.isdigit() else c for c in re.split(r"([0-9]+)", s)]


def sort_nicely(names):
    """Sort fault.00001.nc ... fault.04483.nc numerically, in place."""
    names.sort(key=_alphanum_key)
    return names
