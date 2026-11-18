"""Guards for rule 3a in CI: a skip fails the run unless allow-listed."""
import sys

import pytest

from conftest import ROOT, read

pytestmark = pytest.mark.contract
sys.path.insert(0, str(ROOT / "tests"))
import skip_gate  # noqa: E402

ALLOW = ([("tests/unit/test_lapusta_xi.py", "no gfortran")],
         [("tests/e2e/test_benchmarks.py", "this reference has no summary.json")])


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
    for path, reason in skip_gate.load_allowlist()[0]:
        assert "no bin/eqquasi" not in reason and "no binary" not in reason, \
            (path, reason)
    for path, reason in skip_gate.load_allowlist()[1]:
        assert "binary" not in reason and "gfortran" not in reason, (path, reason)


def test_ci_builds_before_the_fast_tier_and_runs_the_skip_gate():
    wf = read(".github/workflows/test.yml")
    fast = wf.find('-m "not e2e"')
    assert 0 < wf.find("install.eqquasi.sh") < fast, \
        "the fast tier runs before the build, so the binary check cannot pass"
    assert "-rfEs" in wf and "tests/skip_gate.py" in wf, \
        "CI does not print skip reasons and gate on them"


def test_trailing_output_after_the_summary_does_not_disable_the_check():
    line = "SKIPPED [1] tests/contract/test_build_scripts.py:136: no bin/eqquasi-1.24.4"
    log = _log(line) + "\n\n[daemon] prrte exited\n\n"
    assert skip_gate.forbidden_skips(log, ALLOW) == [line]


def test_skip_count_without_reasons_and_trailing_lines_raises():
    with pytest.raises(ValueError, match="-rs"):
        skip_gate.forbidden_skips(_log() + "\n\ndaemon line\n", ALLOW)


def test_summary_text_inside_a_message_is_not_a_summary():
    with pytest.raises(ValueError, match="no pytest summary"):
        skip_gate.forbidden_skips("note: 3 passed, 1 skipped in the plan", ALLOW)


def test_truncated_log_raises():
    with pytest.raises(ValueError, match="no pytest summary"):
        skip_gate.forbidden_skips("SKIPPED [1] tests/a.py:1: x\n3 passed, 1 skipp", ALLOW)


def test_data_absence_skip_is_reported_not_forbidden():
    line = "SKIPPED [2] tests/e2e/test_benchmarks.py:126: this reference has no summary.json"
    log = _log(line, summary="3 passed, 2 skipped in 1.0s")
    assert skip_gate.classify_skips(log, ALLOW) == ([], [line])


def test_data_absence_entry_does_not_cover_another_file():
    line = "SKIPPED [1] tests/unit/test_x.py:3: this reference has no summary.json"
    assert skip_gate.forbidden_skips(_log(line), ALLOW) == [line]


def test_main_exits_zero_for_data_absence_and_one_for_no_binary(tmp_path, capsys):
    absent = "SKIPPED [1] tests/e2e/test_benchmarks.py:97: this reference holds no offfault files"
    nobin = "SKIPPED [1] tests/e2e/test_benchmarks.py:645: no bin/eqquasi-1.24.4"
    ok, bad = tmp_path / "ok.txt", tmp_path / "bad.txt"
    ok.write_text(_log(absent))
    bad.write_text(_log(absent, nobin, summary="3 passed, 2 skipped in 1.0s"))
    assert skip_gate.main([str(ok)]) == 0
    assert skip_gate.main([str(bad)]) == 1
    assert "data-absence skip" in capsys.readouterr().out
