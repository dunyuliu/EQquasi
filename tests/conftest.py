"""Shared fixtures and helpers for the EQquasi regression suite.

These tests are static/structural: they run in seconds and need no MPI, no
MUMPS and no simulation output. They exist so that every defect found by hand
has a mechanical guard, per PROJECT_RULES.md.

Run with:  python3 -m pytest tests/ -v
"""

import pathlib
import re

ROOT = pathlib.Path(__file__).resolve().parent.parent


def read(relpath):
    return (ROOT / relpath).read_text()


def compset_dirs():
    """Every case_input/<name>/ directory holding a user_defined_params.py."""
    base = ROOT / "case_input"
    return sorted(
        d for d in base.iterdir()
        if d.is_dir() and (d / "user_defined_params.py").is_file()
    )


def strip_fortran_comments(text):
    """Drop full-line and trailing Fortran comments, and commented-out blocks."""
    out = []
    for line in text.splitlines():
        stripped = line.lstrip()
        if stripped.startswith("!"):
            continue
        # A '!' outside a quoted string starts a comment.
        depth = None
        for i, ch in enumerate(line):
            if ch in "'\"":
                if depth is None:
                    depth = ch
                elif depth == ch:
                    depth = None
            elif ch == "!" and depth is None:
                line = line[:i]
                break
        out.append(line)
    return "\n".join(out)


def strip_python_comments(text):
    out = []
    for line in text.splitlines():
        if line.lstrip().startswith("#"):
            continue
        out.append(line)
    return "\n".join(out)


def par_attributes_read_by(script_relpath):
    """Every `par.<name>` referenced by a script, comments removed."""
    src = strip_python_comments(read(script_relpath))
    return set(re.findall(r"\bpar\.([A-Za-z_][A-Za-z0-9_]*)", src))


def load_case_params(case):
    """Execute case_input/<case>/user_defined_params.py and return its `par`.

    The case files do `from defaultParameters import parameters`, which lives in
    scripts/, and they build par.on_fault_vars at import time. Executing rather
    than parsing is the point: the guard must see the array the run will
    actually use, including any values computed in the loop.
    """
    import importlib.util
    import sys

    path = ROOT / "case_input" / case / "user_defined_params.py"
    if not path.is_file():
        raise FileNotFoundError(path)

    added = [str(ROOT / "scripts"), str(path.parent)]
    sys.path[:0] = added
    # A previously loaded case would otherwise be returned from the cache.
    stale = sys.modules.pop("user_defined_params", None)
    # Importing must not leave __pycache__/ behind in case_input/<case>/:
    # create.newcase copies that directory's contents file by file, so a stray
    # subdirectory there aborts case creation.
    dont_write = sys.dont_write_bytecode
    sys.dont_write_bytecode = True
    try:
        spec = importlib.util.spec_from_file_location("user_defined_params", path)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        return mod.par
    finally:
        sys.dont_write_bytecode = dont_write
        sys.modules.pop("user_defined_params", None)
        if stale is not None:
            sys.modules["user_defined_params"] = stale
        for p in added:
            if p in sys.path:
                sys.path.remove(p)
