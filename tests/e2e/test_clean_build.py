"""End-to-end: build the solver from clean, the way a new user would.

This is the one thing `pytest -m e2e` covered that nothing else did. Every other
e2e test assumes `bin/eqquasi` already exists and skips when it does not, so on
a fresh checkout they would all skip and the suite would pass without compiling
a line.

`install.eqquasi.sh` is the documented entry point (README, "Build and
installation"), so it is what gets tested -- not `make -C src` directly, which
is what a developer runs and therefore the path least likely to rot.

Opt-in, because a clean build takes minutes:

    python3 -m pytest -m e2e tests/test_clean_build.py
"""

import os
import shutil
import subprocess

import pytest

import cases as C
from conftest import ROOT

pytestmark = [pytest.mark.e2e, pytest.mark.e2e_fast]

INSTALL = ROOT / "install.eqquasi.sh"

BUILD_TIMEOUT_S = 1800


def expected_binary():
    """Where install.eqquasi.sh will put the binary it is about to build.

    bin/ holds versioned binaries and no plain `eqquasi`, so this has to be
    derived from the version src declares rather than hard-coded. It was
    hard-coded to bin/eqquasi, and after the switch to versioned names that
    path stopped existing -- which made the version-skew test below skip
    itself with "no binary" on every run. A guard against stale binaries that
    quietly does not run is worse than no guard, because the suite still
    reports green.
    """
    v = C.declared_version()
    assert v, "EQQUASI_VERSION not found in src/globalvar.f90"
    return ROOT / "bin" / f"eqquasi-{v}"


@pytest.mark.skipif(not INSTALL.exists(), reason="no install.eqquasi.sh")
def test_install_script_produces_a_working_binary():
    """A clean build must yield a binary that runs and reports its version.

    The binary is rebuilt in place, so this deliberately runs before the other
    e2e tests would want it -- pytest orders by file, and a rebuilt binary is
    what they should be testing anyway.
    """
    # Keep the existing binary: if the build fails we restore it rather than
    # leaving the tree without one, which would cascade into every other e2e
    # test skipping and the failure looking like "nothing ran".
    BIN = expected_binary()
    saved = None
    if BIN.exists():
        saved = str(BIN) + ".pretest"
        shutil.copy2(str(BIN), saved)

    env = dict(os.environ)
    env["EQQUASIROOT"] = str(ROOT)
    machine = env.get("MACHINE", "ubuntu")

    try:
        r = subprocess.run(["bash", str(INSTALL), "-m", machine],
                           cwd=str(ROOT), env=env, capture_output=True,
                           text=True, timeout=BUILD_TIMEOUT_S)
        # Name MACHINE in the message: the default is `ubuntu` because that is
        # what CI runs, but on a host whose MUMPS headers sit elsewhere the
        # build dies on a missing dmumps_struc.h, which reads as a code
        # failure rather than as "you did not say which machine this is".
        assert r.returncode == 0, (
            f"install.eqquasi.sh -m {machine} failed (exit {r.returncode}). "
            f"MACHINE was {'set' if 'MACHINE' in os.environ else 'unset, so it defaulted'} "
            f"to {machine}; set MACHINE=<host> if that is the wrong target.\n"
            f"stdout tail:\n{r.stdout[-3000:]}\n"
            f"stderr tail:\n{r.stderr[-3000:]}")
        assert BIN.exists(), f"build reported success but {BIN} is missing"

        # It must actually load -- a binary that cannot find its shared
        # libraries builds fine and fails at run time, which has happened here
        # (a stale EQquasi elsewhere on PATH shadowing this one).
        v = subprocess.run([str(BIN)], cwd=str(ROOT), env=env,
                           capture_output=True, text=True, timeout=120)
        out = v.stdout + v.stderr
        assert "Welcome to EQquasi" in out, (
            "the freshly built binary did not start:\n" + out[-2000:])
    except Exception:
        if saved and os.path.exists(saved):
            shutil.copy2(saved, str(BIN))
        raise
    finally:
        if saved and os.path.exists(saved):
            os.remove(saved)


@pytest.mark.e2e_fast
@pytest.mark.skipif(not INSTALL.exists(), reason="no install.eqquasi.sh")
def test_built_binary_version_matches_the_source():
    """The binary on disk must be built from the current globalvar.f90.

    A stale binary is the quiet way an e2e suite tests last week's code: every
    comparison passes because the oracle was frozen from the same stale build.
    """
    declared = C.declared_version()
    assert declared, "EQQUASI_VERSION not found in src/globalvar.f90"
    BIN = expected_binary()
    assert BIN.exists(), (
        f"src declares {declared} but {BIN.name} is not in bin/. Every e2e "
        f"comparison below would be running some other version. Rebuild: "
        f"EQQUASIROOT=$(pwd) MACHINE=<host> bash install.eqquasi.sh -m <host>")

    r = subprocess.run([str(BIN)], cwd=str(ROOT), capture_output=True,
                       text=True, timeout=120)
    out = r.stdout + r.stderr
    assert declared in out, (
        f"{BIN.name} reports a different version than src declares "
        f"({declared}); the file name is not evidence of what is inside it.")
