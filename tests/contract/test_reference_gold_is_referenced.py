"""Guard for PROJECT_RULES.md rule 8a: every reference/gold file has a reader.

An oracle nobody reads is dead weight that looks like a safety net. This does
not (and cannot) prove a file is compared *correctly* -- only that something in
the test suite or the plotting tooling names it, so an unused file is at least
visible rather than silently accumulating.

Scoped to reference/<bench>/gold/*.{json,nc,csv} -- the oracles rule 3/8 are
about. reference/<bench>/plots/ and reference/bp8/archive/ are byproducts and
scratch, not oracles, and are deliberately excluded.
"""

import sys

from conftest import ROOT, read

GOLD_DIRS = sorted((ROOT / "reference").glob("*/gold"))

SEARCH_ROOTS = [
    "tests",
    "scripts",
    "check.test.py",
    "testNameList.py",
]


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


def _bp8_dynamic_names():
    """BP8's station/profile gold filenames are built at runtime by string
    formatting in tests/e2e/test_bp8_against_gold.py, so no literal substring
    match is possible for them. Import that module's own generators instead of
    guessing -- this is the one place a plain text search cannot see what a
    test actually reads.
    """
    sys.path.insert(0, str(ROOT / "tests" / "e2e"))
    try:
        import test_bp8_against_gold as m
        names = {f"fltst_strk{s}.csv" for s in m.ALL_STATIONS}
        names |= {f.replace(".dat", ".csv") for f in m.PROFILE_FILES}
        return names
    finally:
        sys.path.remove(str(ROOT / "tests" / "e2e"))
        sys.modules.pop("test_bp8_against_gold", None)


def gold_files():
    out = []
    for d in GOLD_DIRS:
        for f in d.glob("*"):
            if f.is_file() and f.suffix in (".json", ".nc", ".csv"):
                out.append(f)
    return out


def test_every_gold_file_is_named_somewhere_in_the_suite():
    hay = _haystack()
    dynamic = _bp8_dynamic_names()
    unreferenced = [
        f.relative_to(ROOT) for f in gold_files()
        if f.name not in hay and f.name not in dynamic
    ]
    assert not unreferenced, (
        "these reference/gold files are not named by any test, script, or "
        "check.test.py/testNameList.py, and do not match a BP8 dynamically "
        "generated name either -- either wire them into a test in the same "
        "change that adds them, or do not commit them (PROJECT_RULES.md rule "
        f"8a): {unreferenced}"
    )
