"""Guards for the regression gate itself.

The gate was vacuous: check.test.py printed `FAIL ` on divergence, but CI's
grep pattern did not match that string and testAll.py swallowed the exit code
with `|| true`. A BP5 regression would have gone green. PROJECT_RULES.md rule 3
names this the highest-priority mechanical gap.
"""

import re
import subprocess
import sys

from conftest import ROOT, read


WORKFLOW = ".github/workflows/test.yml"


def test_ci_pattern_matches_the_failure_string_check_test_actually_prints():
    """The exact string check.test.py emits must be caught by CI's grep."""
    check = read("check.test.py")
    assert "'FAIL '" in check, "check.test.py no longer prints the literal 'FAIL '"

    wf = read(WORKFLOW)
    patterns = re.findall(r"grep -q(?:i)?E ['\"]([^'\"]+)['\"]", wf)
    assert patterns, "no grep pattern found in the workflow"

    sample = "FAIL reference/bp5/fault.00101.nc test/x/Q0/fault.00101.nc"
    matched = any(
        subprocess.run(["grep", "-qE", p], input=sample, text=True).returncode == 0
        or subprocess.run(["grep", "-qiE", p], input=sample, text=True).returncode == 0
        for p in patterns
    )
    assert matched, (
        "no CI grep pattern matches check.test.py's FAIL output; "
        f"patterns were {patterns}"
    )


def test_ci_does_not_swallow_the_test_exit_code():
    wf = read(WORKFLOW)
    assert "testAll.py >> testRunLog.txt || true" not in wf, \
        "`|| true` discards the test suite's verdict"
    assert "|| true" not in wf


def test_ci_refuses_to_pass_if_the_check_never_ran():
    """A crash before check.test.py must not read as success."""
    wf = read(WORKFLOW)
    assert "CHECK SUMMARY" in wf, \
        "CI should require positive evidence that check.test.py completed"


def test_check_test_exits_nonzero_on_failure():
    """check.test.py must signal failure through its exit status, not just stdout."""
    src = read("check.test.py")
    assert "sys.exit" in src, "check.test.py never sets a non-zero exit status"


def test_check_test_does_not_silently_skip_a_missing_output():
    """A reference that exists but an output that does not is a failure, not a skip."""
    src = read("check.test.py")
    assert "missing output" in src, \
        "a produced-file-missing case must be reported as FAIL, not skipped"


def test_check_test_fails_when_nothing_was_compared():
    src = read("check.test.py")
    assert "vacuous" in src or "compared == 0" in src, \
        "a run that compares zero files must not pass"


def test_testall_propagates_the_verdict():
    src = read("testAll.py")
    assert "sys.exit" in src, "testAll.py must exit non-zero when the check fails"


def test_testall_has_no_dead_mpirun_variable():
    """MPIRUN was defined and never used, which misleads anyone debugging launch."""
    src = read("testAll.py")
    assert "MPIRUN" not in src


def test_check_test_gate_is_self_consistent_end_to_end(tmp_path):
    """Run check.test.py where every output is missing: it must fail.

    Done in an isolated tree, not in ROOT, so the verdict does not depend on
    whether a previous run happened to leave a populated test/ directory
    behind. An earlier version of this test read ROOT directly and silently
    flipped to passing once the physics suite had been run.
    """
    import shutil

    for name in ("check.test.py", "testNameList.py"):
        shutil.copy(ROOT / name, tmp_path / name)
    # References moved to reference/<benchmark>/ when the BP5/BP7 and BP8
    # oracles were unified; check.test.py reads them from there.
    shutil.copytree(ROOT / "reference", tmp_path / "reference")
    # Deliberately do NOT create test/ -- references exist, outputs do not.

    r = subprocess.run(
        [sys.executable, "check.test.py"],
        cwd=str(tmp_path), capture_output=True, text=True,
    )
    assert r.returncode != 0, (
        "check.test.py passed with no simulation output present:\n" + r.stdout
    )
    assert "missing output" in r.stdout, (
        "a missing output should be reported as FAIL, not skipped:\n" + r.stdout
    )
