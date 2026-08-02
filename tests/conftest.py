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
