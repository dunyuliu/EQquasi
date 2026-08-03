"""Guards for where run artifacts are allowed to land.

Convention: every scratch artifact -- generated cases, simulation output, build
products -- lives under `work/` at the repo root, which is gitignored. Nothing
scratch is written to the repo root itself.

This matters because `testAll.py` starts by deleting its scratch directory. A
scratch directory at the repo root is one typo away from deleting something
tracked, and it also means `git status` is never clean after a run.
"""

import re

from conftest import ROOT, read

SCRATCH_ROOT = "work"


def gitignore_entries():
    return {ln.strip() for ln in read(".gitignore").splitlines() if ln.strip()}


def test_scratch_root_is_gitignored():
    entries = gitignore_entries()
    assert SCRATCH_ROOT + "/" in entries or SCRATCH_ROOT in entries, \
        f"{SCRATCH_ROOT}/ must be gitignored; run artifacts are never committed"


def test_test_harness_writes_under_the_scratch_root():
    """testAll.py used to create and rm -rf a `test/` directory at the repo root."""
    src = read("testAll.py")
    assert "rm -rf work/test" in src, \
        "testAll.py must clear its scratch dir under work/, not at the repo root"
    assert re.search(r"rm -rf ['\"]?test['\"]?[\s'\"]*\)", src) is None, \
        "testAll.py must not rm -rf a bare `test` directory at the repo root"
    assert "work/test" in src


def test_check_test_reads_from_the_scratch_root():
    src = read("check.test.py")
    m = re.search(r"testRoot\s*=\s*['\"]([^'\"]+)['\"]", src)
    assert m, "could not find testRoot in check.test.py"
    assert m.group(1).startswith(SCRATCH_ROOT + "/"), \
        f"check.test.py must read outputs from {SCRATCH_ROOT}/, got {m.group(1)!r}"


def test_reference_results_are_not_under_the_scratch_root():
    """The oracles are committed data and must survive a scratch wipe."""
    src = read("check.test.py")
    m = re.search(r"refRoot\s*=\s*['\"]([^'\"]+)['\"]", src)
    assert m, "could not find refRoot in check.test.py"
    assert not m.group(1).startswith(SCRATCH_ROOT + "/"), \
        "reference results must NOT live under the gitignored scratch root"
    assert (ROOT / m.group(1)).is_dir()


def test_no_tracked_files_under_the_scratch_root():
    """A tracked file under work/ would be destroyed by the next test run."""
    import subprocess
    r = subprocess.run(["git", "ls-files", SCRATCH_ROOT],
                       cwd=str(ROOT), capture_output=True, text=True)
    tracked = [ln for ln in r.stdout.splitlines() if ln.strip()]
    assert not tracked, (
        f"these files are tracked under {SCRATCH_ROOT}/ and would be deleted by "
        f"the next run: {tracked}"
    )
