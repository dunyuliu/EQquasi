"""Guards for rule 3a in CI: a skip fails the run unless allow-listed."""
import sys

import pytest

from conftest import ROOT, read

pytestmark = pytest.mark.contract
sys.path.insert(0, str(ROOT / "tests"))
import skip_gate  # noqa: E402

ALLOW = [("tests/unit/test_lapusta_xi.py", "no gfortran")]


def _log(*skips, summary="3 passed, 1 skipped in 1.0s"):
    return "\n".join(["...", "=== short test summary info ===", *skips, summary])


def test_a_no_binary_skip_is_forbidden():
    line = "SKIPPED [1] tests/contract/test_build_scripts.py:136: no bin/eqquasi-1.24.4"
    assert skip_gate.forbidden_skips(_log(line), ALLOW) == [line]


def test_an_allow_listed_host_tool_skip_passes():
    line = "SKIPPED [1] tests/unit/test_lapusta_xi.py:61: no gfortran"
    assert skip_gate.forbidden_skips(_log(line), ALLOW) == []


def test_allow_list_needs_file_and_reason_to_match():
    line = "SKIPPED [1] tests/unit/test_other.py:9: no gfortran"
    assert skip_gate.forbidden_skips(_log(line), ALLOW) == [line]


def test_log_without_summary_raises():
    with pytest.raises(ValueError, match="no pytest summary"):
        skip_gate.forbidden_skips("collected 3 items", ALLOW)


def test_skip_count_without_reasons_raises():
    with pytest.raises(ValueError, match="-rs"):
        skip_gate.forbidden_skips(_log(), ALLOW)


def test_malformed_allow_list_line_raises(tmp_path):
    p = tmp_path / "a.txt"
    p.write_text("tests/x.py no tab here\n")
    with pytest.raises(ValueError, match="expected"):
        skip_gate.load_allowlist(p)


def test_committed_allow_list_never_excuses_a_missing_binary():
    for path, reason in skip_gate.load_allowlist():
        assert "no bin/eqquasi" not in reason and "no binary" not in reason, \
            (path, reason)


def test_ci_builds_before_the_fast_tier_and_runs_the_skip_gate():
    wf = read(".github/workflows/test.yml")
    fast = wf.find('-m "not e2e"')
    assert 0 < wf.find("install.eqquasi.sh") < fast, \
        "the fast tier runs before the build, so the binary check cannot pass"
    assert "-rs" in wf and "tests/skip_gate.py" in wf, \
        "CI does not print skip reasons and gate on them"
