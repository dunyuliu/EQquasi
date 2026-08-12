import pytest
"""Guard for PROJECT_RULES.md rule 12: the gate set must include ntotft > 1.

BP5, BP5-dip90, BP7 and BP8 all have ntotft = 1 (case_input/test.bp5.qdc,
test.bp5.qdc.dip90, test.bp7.qdc, test.bp8.qdc). Three multi-fault bugs in this
project's history were found by hand, running a new case -- none were caught by
the gate, because nothing in it ever exercised ntotft > 1. This is expected to
FAIL until a multi-fault compset (case_input/test.stepover.qdc and
case_input/bp1002.qdc.2000 both already have ntotft = 2) is added to
the CASES table in tests/e2e/cases.py, and gets its own
reference/<bench>/.
"""

import sys

from conftest import ROOT, compset_dirs, load_case_params, read


def gate_compset_names():
    """Every compset the regression/e2e gate actually runs."""
    names = set()

    # The gate set is the CASES table in tests/e2e/cases.py; testNameList.py
    # was deleted when the five per-benchmark e2e files were collapsed.
    sys.path.insert(0, str(ROOT / "tests" / "e2e"))
    try:
        import cases as C
        names.update(c[1] for c in C.CASES)
    finally:
        sys.path.remove(str(ROOT / "tests" / "e2e"))

    e2e_dir = ROOT / "tests" / "e2e"
    for f in e2e_dir.glob("*.py"):
        src = read(f"tests/e2e/{f.name}")
        for m in (str(d.relative_to(ROOT / "case_input")) for d in compset_dirs()):
            if f'"{m}"' in src or f"'{m}'" in src:
                names.add(m)
    return names


@pytest.mark.xfail(reason=(
    "The step-over reference was removed in v1.10.0 along with BP5-dip90 and "
    "BP7: produced ad hoc, before the workflow rule. It is to be recreated "
    "through the e2e mechanism, at which point this passes again. Kept as "
    "xfail rather than deleted because the gap is real -- three multi-fault "
    "bugs reached master while every gate ran a single fault."), strict=False)
def test_gate_set_includes_a_multifault_compset():
    import sys
    sys.path.insert(0, str(ROOT))
    try:
        names = gate_compset_names()
    finally:
        sys.path.remove(str(ROOT))

    multifault = []
    for name in sorted(names):
        try:
            par = load_case_params(name)
        except FileNotFoundError:
            continue
        if getattr(par, "ntotft", 1) > 1:
            multifault.append(name)

    assert multifault, (
        "every compset the regression/e2e gate runs "
        f"({sorted(names)}) has ntotft == 1; three multi-fault bugs in this "
        "project were found only by running a new case by hand. Add "
        "test.stepover.qdc or bp1002.qdc.2000 (both ntotft = 2) to "
        "the CASES table in tests/e2e/cases.py, with its own "
        "reference/<bench>/ (PROJECT_RULES.md rule 12)."
    )
