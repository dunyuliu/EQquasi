"""Guard: CI must actually gate on the e2e result.

The gate was vacuous once already. `check.test.py` printed `FAIL ` on
divergence, but CI's grep pattern did not match that string and `testAll.py`
swallowed the exit code, so a run where nothing matched its reference reported
success. Both files are gone now -- `pytest -m e2e_fast` replaced them -- and
the failure mode has not: a CI step that runs tests but does not fail on them
is worse than no step, because it looks like coverage.

So this asserts the three things that make the gate real:

  1. CI runs the fast e2e tier, not just the static tiers.
  2. Nothing suppresses the exit code (`|| true`, a bare `tee` without
     `set -o pipefail`).
  3. The tiers exist and are non-empty, so `-m e2e_fast` cannot silently
     select nothing.
"""

import os
import subprocess
import sys

import pytest

from conftest import ROOT, read

pytestmark = pytest.mark.contract

WORKFLOW = os.path.join(".github", "workflows", "test.yml")


def _wf():
    p = ROOT / WORKFLOW
    if not p.exists():
        pytest.skip("no GitHub workflow")
    return read(WORKFLOW)


def test_ci_runs_the_fast_e2e_tier():
    wf = _wf()
    assert "-m e2e_fast" in wf, (
        "CI does not run the fast e2e tier. The static tiers read source and "
        "never execute the solver, so on their own they cannot catch a change "
        "that alters the physics while leaving the text intact.")


def test_ci_does_not_swallow_the_exit_code():
    wf = _wf()
    bad = [ln.strip() for ln in wf.splitlines()
           if "pytest" in ln and ("|| true" in ln or "|| :" in ln)]
    assert not bad, f"CI suppresses a pytest failure: {bad}"

    # `cmd | tee log` reports tee's status, not cmd's, unless pipefail is set.
    for ln in wf.splitlines():
        if "pytest" in ln and "tee" in ln:
            assert "set -o pipefail" in wf, (
                "a pytest command is piped into tee without `set -o pipefail`, "
                "so the pipeline reports tee's exit status and a failing suite "
                "passes CI")


def test_both_e2e_tiers_select_something():
    """`-m e2e_fast` selecting zero tests would pass CI while testing nothing."""
    for marker in ("e2e_fast", "e2e"):
        r = subprocess.run(
            [sys.executable, "-m", "pytest", "tests/", "-q", "-m", marker,
             "--collect-only"],
            cwd=str(ROOT), capture_output=True, text=True, timeout=300)
        import re as _re
        m = _re.search(r"(\d+)\s*/\s*\d+ tests collected", r.stdout) or \
            _re.search(r"(\d+) tests? collected", r.stdout)
        assert m, f"could not read a collection count from:\n{r.stdout[-500:]}"
        assert int(m.group(1)) > 0, \
            f"-m {marker} selects no tests; the gate would be vacuous"


def test_the_case_table_covers_both_tiers():
    """A tier with no case in the table is a tier that tests nothing."""
    sys.path.insert(0, str(ROOT / "tests" / "e2e"))
    try:
        import cases as C
    finally:
        sys.path.remove(str(ROOT / "tests" / "e2e"))
    tiers = {c[4] for c in C.CASES}
    assert "fast" in tiers, "no fast case in tests/e2e/cases.py"
    assert "full" in tiers, "no full-only case in tests/e2e/cases.py"
