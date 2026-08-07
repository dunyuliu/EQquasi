"""End-to-end: BP5 and BP7 against the committed reference results.

This tier already existed for BP5 and BP7, as `testAll.py` plus `check.test.py`
diffing `fault.00101.nc` against `test.reference.results/`. It simply lived
outside pytest, so the suite's tiers were split across two runners and
`python3 -m pytest tests/` silently covered no physics at all.

This wraps it so all three benchmarks report through one command. It does not
reimplement the comparison -- `check.test.py` remains the single definition of
what "matches the reference" means, and its exit code is the verdict.

Note what the references are: *regression locks*, not physics validations. The
test compsets are deliberately coarse (`test.bp5.qdc` uses dx = 4000 m against
BP5's 2000 m; `test.bp7.qdc` uses dx = 50 m, about four cells across BP7's
velocity-weakening disc). They catch unintended change. They do not establish
that a benchmark is reproduced correctly -- that needs the production compset at
its native resolution, and for BP8 it needs `reference/bp8/gold/`.

Opt-in with the rest of the e2e tier:

    python3 -m pytest -m e2e tests/

`testAll.py` recompiles from source and wipes `work/test/` before running, so do
not launch it alongside another build.
"""

import os
import subprocess

import pytest

from conftest import ROOT

RUN_TIMEOUT_S = 7200


@pytest.mark.e2e
def test_bp5_and_bp7_reproduce_reference_results():
    """testAll.py builds, runs three compsets and diffs against the oracles."""
    refs = os.path.join(str(ROOT), "test.reference.results")
    if not os.path.isdir(refs):
        pytest.skip("no test.reference.results/ to compare against")

    env = dict(os.environ)
    env["EQQUASIROOT"] = str(ROOT)
    env["PATH"] = f"{ROOT}/bin:{ROOT}/scripts:" + env.get("PATH", "")
    env["OMP_NUM_THREADS"] = "1"

    try:
        r = subprocess.run(["python3", "testAll.py"], cwd=str(ROOT), env=env,
                           capture_output=True, text=True, timeout=RUN_TIMEOUT_S)
    except subprocess.TimeoutExpired:
        pytest.fail(f"testAll.py exceeded {RUN_TIMEOUT_S} s")

    # check.test.py exits non-zero on a regression and testAll.py propagates it.
    assert r.returncode == 0, (
        "BP5/BP7 regression against test.reference.results/ failed.\n"
        f"--- stdout tail ---\n{r.stdout[-4000:]}\n"
        f"--- stderr tail ---\n{r.stderr[-2000:]}"
    )
    assert "FAIL" not in r.stdout, (
        f"check.test.py printed FAIL while exiting 0:\n{r.stdout[-4000:]}"
    )
