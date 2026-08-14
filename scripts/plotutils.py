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


# Figure style lives once, at the bottom of this file under "aesthetics", and
# is applied on import. There used to be a second apply_style() up here: the
# module-level call ran THIS one, the later def rebound the name, and every
# utility's explicit pu.apply_style() then partially overrode it -- so the
# style actually in force was a merge of the two that neither one stated.


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
    # Cycles live under result/ since the case gained input/, result/,
    # scratch/, tool/ and log/. Older cases kept them at case level, and older
    # still called them Q0, Q1 -- accept all three so existing output stays
    # readable rather than needing every case regenerated.
    for d in glob.glob(os.path.join(path, "result", "cycle[0-9]*")) + \
             glob.glob(os.path.join(path, "cycle[0-9]*")) + \
             glob.glob(os.path.join(path, "Q[0-9]*")):
        m = re.fullmatch(r"(?:cycle|Q)(\d+)", os.path.basename(d))
        if m and os.path.isdir(d) and ok(d):
            qs.append((int(m.group(1)), d))
    if qs:
        return [(os.path.basename(d), d) for _, d in sorted(qs)]
    die(f"{tag}found none of the required files ({', '.join(patterns)}) in "
        f"{path}, in {os.path.join(path, 'result', 'cycle*')}, or in any "
        f"cycle*/Q* directory beside it.\n"
        "Point me at a cycle directory, or at a case directory whose "
        "result/cycle* folders hold the run output.")


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
            hit = next((c for c in
                        (f"result/cycle{int(m.group(1))}" if m else "",
                         f"cycle{int(m.group(1))}" if m else "",
                         f"Q{int(m.group(1))}" if m else "")
                        if c and os.path.isdir(c)), None)
            if hit:
                cand = hit
            else:
                have = sorted(os.path.basename(d)
                              for d in glob.glob("result/cycle[0-9]*")
                              + glob.glob("cycle[0-9]*") + glob.glob("Q[0-9]*")
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
    """The `par` object from user_defined_params.py at or above run_dir.

    A cycle folder holds output only; the parameters live in the case root,
    which is two levels up now that cycles sit in result/cycleN (one level for
    older cases with cycleN at case level). Loads by explicit path -- never via
    cwd or PYTHONPATH -- with the case's tool/ and scripts/ on sys.path for
    `from defaultParameters import ...`.
    """
    import importlib.util
    tag = f"{tool}: " if tool else ""
    run_dir = os.path.abspath(run_dir or ".")
    up1 = os.path.dirname(run_dir)
    for d in (run_dir, up1, os.path.dirname(up1)):
        p = os.path.join(d, "user_defined_params.py")
        if os.path.isfile(p):
            break
    else:
        die(f"{tag}no user_defined_params.py in {run_dir} or the two "
            "directories above it; this tool needs the case parameters "
            "(point it at the case or cycle directory).")
    added = [os.path.join(d, "tool"),
             os.path.dirname(os.path.abspath(__file__)),
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


def _depth_label(par, zlo):
    """"z (km)", saying how far below the fault the domain goes.

    A fault-plane map is drawn to the fault's own extent, so a fault that
    bottoms out at -20 km looks identical whether the domain below it is
    20 km deep or 60 km. That matters here: the depth of the domain under
    the seismogenic zone is what the loading is applied through, and reading
    these panels as if the model ended at the fault's base is a mistake the
    axis should not invite.

    Falls back to a plain label when the case declares no domain depth, or
    when the fault reaches the bottom of it and there is nothing to say.
    """
    zdom = getattr(par, "fzmin", None)
    if zdom is None or zlo is None or zdom >= zlo - 1.0:
        return "z (km)"
    return f"z (km)   [domain to {zdom / 1e3:.0f} km]"


def fault_coords(par, ndip, nstrike, ift=0):
    """(x, z, xlabel, ylabel) for fault ift's plane, in km.

    A multi-fault case's par.fx spans the union of the segments, not any one
    fault, so it cannot label a single fault's axes -- falling back to node
    indices there loses the very thing a step-over plot is for, namely where
    each segment actually sits along strike. par.faultgeom carries the
    per-fault (xlo, xhi, ycoor, zlo, zhi), so use it when present and only
    then fall back."""
    geom = getattr(par, "faultgeom", None)
    if geom is not None and ift < len(geom):
        xlo, xhi, _ycoor, zlo, zhi = geom[ift]
        return (np.linspace(xlo, xhi, nstrike) / 1e3,
                np.linspace(zlo, zhi, ndip) / 1e3,
                "along strike (km)", _depth_label(par, zlo))
    try:
        fx, fz = np.asarray(par.fx), np.asarray(par.fz)
        if fx.size == nstrike and fz.size == ndip:
            return fx / 1e3, fz / 1e3, "along strike (km)", "z (km)"
    except AttributeError:
        pass
    return (np.arange(nstrike), np.arange(ndip),
            "node # along strike", "node # along dip")


def fault_y(par, ift):
    """Fault ift's fault-normal offset in km, or None if the case does not
    declare one (single-fault cases, legacy compsets)."""
    geom = getattr(par, "faultgeom", None)
    if geom is not None and ift < len(geom):
        return geom[ift][2] / 1e3
    return None


# ----------------------------------------------------------------- aesthetics

def apply_style():
    """The one figure style. Values are what was already in force.

    These are not a redesign: they are the merge the two former definitions
    produced between them, read off matplotlib.rcParams and written down, so
    that consolidating them changed no figure. Font sizes suit a canvas at
    roughly print size; a script that renders larger than it prints must scale
    them by k = canvas width / print width itself, because a point size is
    only meaningful against the width the figure is reproduced at.
    """
    import matplotlib
    matplotlib.use("Agg")
    matplotlib.rcParams.update({
        "figure.dpi": 150,
        "savefig.dpi": 150,
        "font.size": 10,
        "axes.titlesize": 11,
        "axes.labelsize": 10,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "legend.fontsize": 9,
        "figure.titlesize": 13,
        "axes.linewidth": 1.1,
        "lines.linewidth": 2.2,
        "lines.markersize": 8,
        "xtick.major.width": 1.1,
        "ytick.major.width": 1.1,
        "xtick.major.size": 5.5,
        "ytick.major.size": 5.5,
        "grid.alpha": 0.25,
        "axes.grid": False,
        "savefig.bbox": "tight",
    })


apply_style()


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
