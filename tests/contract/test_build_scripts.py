"""Guards for the build and install path.

Every assertion here corresponds to a defect that was hit for real while
following the README to run BP5 on a host without system MUMPS.
"""

import re
import shutil
import subprocess
import sys

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
    """mpif90.openmpi is broken on some hosts; the makefile must not hardcode it.

    Whether a working MPI Fortran compiler exists is a property of the host, not
    of this repository, so skip where none is installed rather than reporting a
    repo defect.
    """
    import pytest
    if not (shutil.which("mpif90") or shutil.which("mpif90.openmpi")):
        pytest.skip("no MPI Fortran wrapper on this host")
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


def test_binary_on_disk_is_built_from_the_current_source():
    """bin/eqquasi's version must match what src/globalvar.f90 declares.

    This is the fast-tier copy of the same check e2e makes before running a
    benchmark, and it is here because the e2e one only speaks after a 14-minute
    run. On 2026-08-12 a whole session of BP1002 step-over work -- runs,
    numbers quoted in a deck, and a part-generated reference -- was produced
    against a 1.7.0 binary while src declared 1.10.0, across commits that
    changed the multi-fault fault-node engine that case exercises. The e2e
    tier had been failing on exactly this and nobody read it in time.

    A stale binary is silent by construction: comparisons pass, because the
    oracle was frozen from the same stale build. Catch it in seconds, on the
    tier that runs on every push.
    """
    import pytest
    binary = ROOT / "bin" / "eqquasi"
    if not binary.exists():
        pytest.skip("no bin/eqquasi; nothing built yet")

    # Share the e2e tier's implementation rather than paraphrasing it. The two
    # copies had drifted into different semantics -- one matched the version as
    # a substring, the other tokenised it; one failed on an unreadable binary,
    # the other passed silently -- so "the check passes" meant different things
    # depending on which tier you asked.
    sys.path.insert(0, str(ROOT / "tests" / "e2e"))
    from cases import binary_version, declared_version

    declared = declared_version()
    assert declared, "EQQUASI_VERSION not found in src/globalvar.f90"
    try:
        built = binary_version(str(binary))
    except Exception as exc:
        pytest.fail(f"could not read bin/eqquasi's version ({exc}); a binary "
                    f"of unknown provenance cannot be trusted to have been "
                    f"built from this source.")
    assert built == declared, (
        f"bin/eqquasi reports {built}, but src/globalvar.f90 declares "
        f"{declared}. Every result and reference produced from it describes "
        f"older code. Rebuild: EQQUASIROOT=$(pwd) MACHINE=<host> make -C src "
        f"&& mv src/eqquasi bin/   (MACHINE=utig on the utig hosts)")
