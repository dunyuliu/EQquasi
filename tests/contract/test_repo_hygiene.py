"""Guards for where run artifacts are allowed to land.

Convention: every scratch artifact -- generated cases, simulation output, build
products -- lives under `work/` at the repo root, which is gitignored. Nothing
scratch is written to the repo root itself.

This matters because `pytest -m e2e` starts by deleting its scratch directory. A
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


def test_e2e_runs_stay_under_the_scratch_root():
    """Runs go in work/, which is gitignored -- never beside the sources.

    Replaces three checks that policed testAll.py and check.test.py. Both are
    gone; cases.run_case is now the single place a run directory is chosen, so
    this asserts the property there instead of in two harnesses.
    """
    bench = read("tests/e2e/test_benchmarks.py")
    m = re.search(r'WORK_ROOT\s*=\s*os\.path\.join\(str\(ROOT\),\s*"([^"]+)"', bench)
    assert m, "test_benchmarks.py no longer names a work root"
    assert m.group(1) == SCRATCH_ROOT, \
        "e2e runs are written outside " + SCRATCH_ROOT + "/"
    assert "shutil.rmtree(workdir" in read("tests/e2e/cases.py"), \
        "cases.run_case no longer clears its work directory, so a stale result "\
        "could be mistaken for a fresh one"


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


# --- version provenance ---------------------------------------------------

def declared_version():
    m = re.search(r"EQQUASI_VERSION\s*=\s*'([^']+)'", read("src/globalvar.f90"))
    assert m, "src/globalvar.f90 must declare EQQUASI_VERSION"
    return m.group(1)


def test_version_has_a_single_source_of_truth():
    """The version was hardcoded in five places and stamped into submissions."""
    import glob
    from conftest import ROOT, strip_fortran_comments

    v = declared_version()
    offenders = []
    for path in glob.glob(str(ROOT / "src" / "*.f90")):
        rel = "src/" + path.rsplit("/", 1)[1]
        src = strip_fortran_comments(read(rel))
        for n, line in enumerate(src.splitlines(), 1):
            if "EQQUASI_VERSION" in line:
                continue
            if re.search(r"['\"][^'\"]*" + re.escape(v) + r"[^'\"]*['\"]", line):
                offenders.append(f"{rel}:{n}")
    assert not offenders, (
        f"version {v!r} is hardcoded at {offenders}; use EQQUASI_VERSION so the "
        "string stamped into benchmark submissions cannot drift"
    )


def test_declared_version_matches_a_git_tag():
    """The code claimed 1.3.3 while the newest tag was v1.3.2.

    That version is written into every benchmark output header and into
    runInfo.json, so it is the provenance for a published comparison. A version
    that names no reachable commit cannot be checked out by a reviewer.
    """
    import subprocess
    from conftest import ROOT

    v = declared_version()
    if v.endswith("-dev"):
        # Unreleased work in progress. The '-dev' suffix is deliberate: it is
        # written into benchmark submission headers, where it truthfully signals
        # to a reviewer that the producing code is not a tagged release.
        return
    r = subprocess.run(["git", "tag"], cwd=str(ROOT), capture_output=True, text=True)
    tags = {t.strip().lstrip("v").replace("_", ".") for t in r.stdout.splitlines()}
    if not tags:
        # A shallow clone (actions/checkout without fetch-depth: 0) has no tags,
        # so there is nothing to check against. Skip rather than claim a defect.
        import pytest
        pytest.skip("no git tags available; shallow clone?")
    assert v in tags, (
        f"EQQUASI_VERSION is {v!r} but no matching git tag exists (tags: "
        f"{sorted(tags)}). Either tag this release, bump the declared version, "
        f"or mark it {v}-dev while it is unreleased."
    )


def test_generated_launchers_pin_openmp_threads():
    """Unset OMP_NUM_THREADS lets each rank spawn unbounded threads.

    The build uses -fopenmp. With the variable unset, several MPI ranks on one
    node oversubscribe it badly -- a load average of 62 on 64 cores was observed
    from eight ranks. Measured MPI scaling is only 1.24x for 8 ranks, so threads
    buy nothing at this problem size and the launcher should pin them.
    """
    src = read("scripts/case.setup")
    launchers = [b for b in src.split("def create_") if "f.write" in b]
    pinning = [b for b in launchers if "OMP_NUM_THREADS" in b]
    assert len(pinning) >= 3, (
        "run.sh and the batch scripts should all export OMP_NUM_THREADS; "
        f"only {len(pinning)} launcher(s) do"
    )
    assert "${OMP_NUM_THREADS:-1}" in src, \
        "pin to 1 by default but let the caller override"


def test_no_dead_aztec_sources_in_src():
    """The Aztec solver path is unreachable; its sources belong in archive/.

    src/eqquasi.f90 prints "aztec is temporarily disabled" for sol_op==2 and the
    call is commented out. The makefile's build rules for main_aztec.o and
    elemcal_aztec.o were commented out too, so those files were never compiled --
    but msr.f90 stayed in OBJ and was compiled and linked on every build despite
    its only caller being main_aztec.f90. All three now live under archive/aztec/.

    AZTEC_OPTIONS is deliberately NOT removed: read_input.f90 still reads it from
    model.txt, which is a positional contract with scripts/case.setup. Dropping
    the field would shift every value after it.
    """
    import os
    from conftest import ROOT

    for name in ("main_aztec.f90", "elemcal_aztec.f90", "msr.f90"):
        assert not os.path.exists(os.path.join(str(ROOT), "src", name)), (
            f"src/{name} is dead code; it belongs in archive/aztec/"
        )

    mk = read("src/makefile")
    for obj in ("msr.o", "main_aztec.o", "elemcal_aztec.o"):
        assert obj not in strip_makefile_comments(mk), (
            f"{obj} is still built, but its source is archived"
        )


def strip_makefile_comments(text):
    return "\n".join(l for l in text.splitlines() if not l.lstrip().startswith("#"))
