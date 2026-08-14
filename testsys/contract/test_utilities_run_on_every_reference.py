"""Guard for PROJECT_RULES.md rule 18: post-processing works on every reference.

Every post-processing utility must run against every frozen reference result and
either produce a figure or explain itself. Neither a traceback nor a silent
no-op is acceptable.

This exists because a whole class of defect only appears when a tool meets a
benchmark it was not written against, and the per-benchmark tests could not see
any of it:

  - `plotOnFaultVars` read `global.dat` with `np.loadtxt`, which chokes on the
    section 4.2 field-name line BP8 writes and BP5 does not. Four utilities had
    the same defect.
  - `plotStations.py` was written against BP8's 11-column station layout and
    mislabelled the 9-column BP5/BP7 one throughout -- slip rate is linear
    there, not log10, and column 7 is effective normal stress, not pore
    pressure. The figure looked entirely plausible.
  - `plotOnFaultInitals` had never run at all (NameError).
  - BP5 writes stations as `.txt`, BP8 as `.dat`; a single-extension glob finds
    nothing on the other and reports it as "no stations".

The test is cheap because it runs the utilities against `reference/`, which is
already in the repo -- no solver, no build.
"""

import os
import subprocess
import sys

import pytest

from conftest import ROOT

pytestmark = pytest.mark.contract

SCRIPTS = ROOT / "script"

# Utilities that take a results directory and are expected to work anywhere.
# Campaign-specific tools (plotDomainSweep, plotICReadings) are excluded by
# design: they read a fixed set of sweep directories, not an arbitrary run.
UTILITIES = ["plotRuptureTime.py", "plotStations.py", "plotPeakSliprateTime.py"]


def reference_results():
    """Every frozen results directory under reference/, as (id, path)."""
    out = []
    for bench in sorted(os.listdir(ROOT / "reference")):
        gold = ROOT / "reference" / bench
        if not gold.is_dir():
            continue
        # gold/ is a results directory itself when it holds station files, and
        # a container when the results live in cycle0/ or Q0/.
        for cand in [gold] + sorted(p for p in gold.iterdir()
                                    if p.is_dir() and p.name != "plots"):
            if any(f.name.startswith(("fltst_strk", "global"))
                   for f in cand.iterdir() if f.is_file()):
                out.append((f"{bench}/{cand.name}"
                            if cand != gold else bench, str(cand)))
    return out


REFS = reference_results()


@pytest.mark.skipif(not REFS, reason="no reference results in the tree")
@pytest.mark.parametrize("util", UTILITIES)
@pytest.mark.parametrize("ref", [r[1] for r in REFS], ids=[r[0] for r in REFS])
def test_utility_runs_on_reference(util, ref, tmp_path):
    r = subprocess.run(
        [sys.executable, str(SCRIPTS / util), ref, "-o", str(tmp_path)],
        capture_output=True, text=True, timeout=600)

    assert r.returncode == 0, (
        f"{util} failed on {ref} (exit {r.returncode}).\n"
        f"stdout:\n{r.stdout[-2000:]}\nstderr:\n{r.stderr[-2000:]}")

    # A traceback with a zero exit status is still a failure -- it means the
    # tool swallowed the error and carried on.
    assert "Traceback" not in r.stderr, (
        f"{util} raised on {ref} without failing:\n{r.stderr[-2000:]}")

    # Either it wrote something, or it said why not. Silence is the failure
    # mode this test exists to catch.
    wrote = list(tmp_path.glob("*.png")) + list(tmp_path.glob("*.gif"))
    assert wrote or r.stdout.strip(), (
        f"{util} on {ref} produced no figure and printed nothing. If there is "
        "genuinely nothing to plot, say so on stdout.")
