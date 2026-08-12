"""Guard for PROJECT_RULES.md rule 8a: every reference/gold file has a reader.

An oracle nobody reads is dead weight that looks like a safety net. This does
not (and cannot) prove a file is compared *correctly* -- only that something in
the test suite or the plotting tooling names it, so an unused file is at least
visible rather than silently accumulating.

Scoped to reference/<bench>/*.{json,nc,csv} -- the oracles rule 3/8 are
about. reference/<bench>/plots/ and reference/bp8/archive/ are byproducts and
scratch, not oracles, and are deliberately excluded.
"""

import sys

from conftest import ROOT, read

# The gold/ layer was removed in v1.9.0: reference/<bench>/ is itself the
# reference, and a run's results may sit in a subdirectory (cycle0,
# cycle0-step101-fast). Take both levels.
GOLD_DIRS = sorted(d for d in (ROOT / "reference").glob("*") if d.is_dir()) + \
            sorted(d for d in (ROOT / "reference").glob("*/*")
                   if d.is_dir() and d.name != "plots" and d.name != "archive")

SEARCH_ROOTS = [
    "tests",
    "scripts",
]

# Reference files read by glob rather than by name. tests/e2e/test_benchmarks.py
# discovers what to compare from what each reference directory contains -- that
# is the point of it, so a new benchmark is a table row rather than a file --
# and a text search cannot see filenames that are never written down. These
# patterns record which families are covered that way.
GLOB_READ = ("fault.", "fltst_strk", "srfst_strk", "global",
             "_strike.", "_depth.", "cplot_", "tdyna", "runInfo")


def _haystack():
    chunks = []
    for root in SEARCH_ROOTS:
        p = ROOT / root
        if p.is_file():
            chunks.append(p.read_text())
        elif p.is_dir():
            for f in p.rglob("*.py"):
                if "__pycache__" in f.parts:
                    continue
                chunks.append(f.read_text())
    return "\n".join(chunks)


def _covered_by_glob(name):
    """True when tests/e2e/test_benchmarks.py would find this file by pattern."""
    return any(name.startswith(pfx) or pfx in name for pfx in GLOB_READ)


def gold_files():
    out = []
    for d in GOLD_DIRS:
        for f in d.glob("*"):
            if f.is_file() and f.suffix in (".json", ".nc", ".csv"):
                out.append(f)
    return out


def test_every_gold_file_is_named_somewhere_in_the_suite():
    hay = _haystack()
    dynamic = set()
    unreferenced = [
        f.relative_to(ROOT) for f in gold_files()
        if f.name not in hay and f.name not in dynamic
        and not _covered_by_glob(f.name)
    ]
    assert not unreferenced, (
        "these reference/gold files are not named by any test, script, or "
        "tests/e2e/test_benchmarks.py, and are not among the families "
        "tests/e2e/test_benchmarks.py discovers by glob -- either wire them "
        "into a test in the same "
        "change that adds them, or do not commit them (PROJECT_RULES.md rule "
        f"8a): {unreferenced}"
    )
