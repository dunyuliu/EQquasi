"""Guards for PROJECT_RULES.md rule 13: known Fortran/shell landmines here.

Each check pins a convention this project adopted the hard way. All three pass
today (verified 2026-08-12); they exist to catch a *regression* back into the
pattern, not to report a currently-open defect.
"""

import glob
import re

from conftest import ROOT, read, strip_fortran_comments

SRC_FILES = sorted(
    "src/" + p.rsplit("/", 1)[1] for p in glob.glob(str(ROOT / "src" / "*.f90"))
)


def test_no_ntotft_equals_one_branching():
    """The fault-node engine must stay neutral to ntotft (git e5364d1, b6010e2).

    `if (ntotft == 1) ... else ...` silently leaves the ntotft > 1 path
    unexercised by every single-fault oracle (BP5/BP5-dip90/BP7/BP8 all have
    ntotft = 1) -- exactly the gap rule 12 is about closing from the other
    direction. One code path, ntotft = 1 the degenerate case.
    """
    pattern = re.compile(r"ntotft\s*(==|\.eq\.)\s*1\b", re.IGNORECASE)
    offenders = []
    for rel in SRC_FILES:
        src = strip_fortran_comments(read(rel))
        for n, line in enumerate(src.splitlines(), 1):
            if pattern.search(line):
                offenders.append(f"{rel}:{n}: {line.strip()}")
    assert not offenders, (
        "ntotft == 1 branch found; the fault-node engine must be neutral to "
        f"ntotft (PROJECT_RULES.md rule 13a): {offenders}"
    )


def test_nint_is_never_called_as_the_intrinsic():
    """globalvar.f90 declares `integer, parameter :: nint = 8`, shadowing the
    intrinsic. Calling nint(x) elsewhere gives "Unclassifiable statement", not
    a type error -- easy to lose an hour to. Use int(... + 0.5d0) instead.
    """
    decl = re.compile(r"\bnint\s*=\s*8\b")
    call = re.compile(r"(?<![A-Za-z0-9_])nint\s*\(")
    offenders = []
    for rel in SRC_FILES:
        src = strip_fortran_comments(read(rel))
        for n, line in enumerate(src.splitlines(), 1):
            if decl.search(line):
                continue  # the parameter declaration itself
            if call.search(line):
                offenders.append(f"{rel}:{n}: {line.strip()}")
    assert not offenders, (
        "nint(...) called as the intrinsic, but nint is a shadowed integer "
        f"parameter (=8) in globalvar.f90; use int(x + 0.5d0) "
        f"(PROJECT_RULES.md rule 13b): {offenders}"
    )


def test_no_output_file_opened_on_unit_6():
    """Unit 6 is stdout. Opening a data file on it silently interleaves
    binary/text output with every WRITE(*,*) diagnostic (or clobbers it).
    """
    pattern = re.compile(r"open\s*\(\s*(unit\s*=\s*)?6\s*,", re.IGNORECASE)
    offenders = []
    for rel in SRC_FILES:
        src = strip_fortran_comments(read(rel))
        for n, line in enumerate(src.splitlines(), 1):
            if pattern.search(line):
                offenders.append(f"{rel}:{n}: {line.strip()}")
    assert not offenders, (
        f"a file is opened on unit 6, which is stdout (PROJECT_RULES.md rule "
        f"13c): {offenders}"
    )
