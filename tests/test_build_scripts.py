"""Guards for the build and install path.

Every assertion here corresponds to a defect that was hit for real while
following the README to run BP5 on a host without system MUMPS.
"""

import re
import shutil
import subprocess

from conftest import ROOT, read


def test_install_script_parses():
    """install.eqquasi.sh must at least be syntactically valid bash."""
    r = subprocess.run(["bash", "-n", str(ROOT / "install.eqquasi.sh")],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr


def test_install_script_local_branch_is_reachable():
    """`elif [ $MACHINE == "local"]` (no space before ]) made the branch dead.

    bash -n does not catch it: it parses, then fails at runtime because the
    test operand becomes the literal string `local]`. Guard the shape instead.
    """
    src = read("install.eqquasi.sh")
    assert '"local"]' not in src, \
        'missing space before "]" in the local-machine test: the branch never fires'
    assert "[!" not in src, \
        'missing space after "[": the file-existence test never fires'


def test_install_script_local_branch_actually_selects():
    """Drive the real conditional shape and confirm 'local' is selected."""
    src = read("install.eqquasi.sh")
    m = re.search(r'elif \[(.*?)\]; then\s*\n\s*MUMPS_LIB_DIR=', src)
    assert m, "could not locate the local-machine branch condition"
    cond = m.group(1)
    r = subprocess.run(
        ["bash", "-c", f'MACHINE=local; if [{cond}]; then echo HIT; else echo MISS; fi'],
        capture_output=True, text=True)
    assert r.stdout.strip() == "HIT", \
        f"local branch not selected for MACHINE=local: {r.stdout!r} {r.stderr!r}"


def test_install_script_mkdir_bin_is_idempotent():
    """`mkdir bin` failed on every rerun after the first install."""
    src = read("install.eqquasi.sh")
    assert re.search(r"\bmkdir\s+bin\b(?!\S)", src) is None or "mkdir -p bin" in src, \
        "use `mkdir -p bin` so reinstalling does not fail"
    assert "mkdir -p bin" in src


def test_install_script_reports_build_failure():
    """The install script had no `set -e`, so a failed build was silent."""
    src = read("install.eqquasi.sh")
    assert ("make ||" in src) or ("set -e" in src), \
        "a failing `make` must abort or report, not fall through to `mv`"


def test_makefile_local_scalapack_comes_from_local_mumps():
    """install.mumps.sh builds ScaLAPACK/BLACS into mumps/build/local/lib.

    The `local` target pointed at the system -lscalapack-openmpi instead, which
    by definition is absent on a host that needed a local MUMPS build.
    """
    src = "\n".join(ln for ln in read("src/makefile").splitlines()
                     if not ln.lstrip().startswith("#"))
    local_block = src.split("else ifeq ($(SYSTEM), local)")[1].split("else ifeq")[0]
    assert "-lscalapack-openmpi" not in local_block, \
        "the local target must not depend on a system scalapack"
    assert "-lscalapack" in local_block and "-lblacs" in local_block


def test_makefile_mpi_wrapper_resolves_to_a_real_executable():
    """mpif90.openmpi is broken on some hosts; the makefile must not hardcode it."""
    src = read("src/makefile")
    assert "MPIFC" in src, "expected a probed MPI compiler wrapper"
    m = re.search(r"MPIFC\s*:=\s*\$\(shell (.*)\)\n", src)
    assert m, "could not find the MPIFC probe"
    r = subprocess.run(["bash", "-c", m.group(1)], capture_output=True, text=True)
    picked = r.stdout.strip()
    assert picked, "MPI wrapper probe produced nothing"
    assert shutil.which(picked), f"probe picked {picked!r}, which is not on PATH"
    # And it must actually work, not merely exist.
    r2 = subprocess.run([picked, "-show"], capture_output=True, text=True)
    assert r2.returncode == 0, f"{picked} exists but fails to run: {r2.stderr}"
