"""Both rate-and-state evolution laws must integrate over the step actually taken.

EQquasi carries two different time steps:

    dt     the fixed CFL minimum from read_input.f90, ~4e-3 s at dx = 50 m
    dtev1  the adaptive quasi-dynamic step, 500 s in the BP8 cases

The state evolution must advance over dtev1. Using dt instead does not fail
loudly -- it silently advances the state by a factor dtev1/dt too little, which
at BP8 settings is ~1e5. Through v1.4.6 rate_state_slip_law did exactly that:

    psi = psiss + (sta - psiss) * dexp(-V2*dt/L)      ! wrong

so psi stayed pinned at its initial value for a whole 23-day run while the
aging law beside it evolved correctly, and the fault locked at V ~ 1e-18 m/s.
The run still completed and still produced plausible-looking slip, which is
what made it hard to spot.

The guard is textual because the failure is textual: the two laws sit in the
same file and must use the same step symbol.
"""

import re

import pytest

from conftest import read, strip_fortran_comments

FRIC = "src/fric.f90"

# Subroutines that integrate a state variable forward one step.
EVOLVING_LAWS = ["rate_state_ageing_law", "rate_state_slip_law"]


def subroutine_body(text, name):
    """Return the source of SUBROUTINE <name>, comments stripped."""
    src = strip_fortran_comments(text)
    m = re.search(
        r"^\s*SUBROUTINE\s+" + name + r"\b(.*?)^\s*end\s+SUBROUTINE\s+" + name + r"\b",
        src,
        re.IGNORECASE | re.DOTALL | re.MULTILINE,
    )
    assert m, f"could not locate SUBROUTINE {name} in {FRIC}"
    return m.group(1)


@pytest.mark.parametrize("law", EVOLVING_LAWS)
def test_state_evolution_uses_adaptive_step(law):
    body = subroutine_body(read(FRIC), law)

    # The exponential relaxation factor is where the step enters.
    decays = re.findall(r"dexp\(\s*-\s*V2\s*\*\s*(dtev1|dt)\s*/\s*L\s*\)", body)
    assert decays, (
        f"{law}: no exp(-V2*<step>/L) relaxation term found. If the integration "
        f"was rewritten, update this guard deliberately rather than deleting it."
    )

    bad = [d for d in decays if d != "dtev1"]
    assert not bad, (
        f"{law} advances the state over '{bad[0]}' instead of 'dtev1'.\n"
        f"dt is the fixed CFL minimum (~4e-3 s); dtev1 is the adaptive step "
        f"(500 s in the BP8 cases). Using dt under-advances the state by ~1e5 "
        f"and silently freezes it -- the run completes and looks plausible."
    )


def test_both_laws_use_the_same_step_symbol():
    """Cross-check: whatever symbol is chosen, the two laws must agree."""
    text = read(FRIC)
    steps = set()
    for law in EVOLVING_LAWS:
        body = subroutine_body(text, law)
        steps.update(re.findall(r"dexp\(\s*-\s*V2\s*\*\s*(dtev1|dt)\s*/\s*L\s*\)", body))
    assert len(steps) == 1, (
        f"the aging and slip laws integrate over different time steps: {sorted(steps)}. "
        f"They advance the same state over the same interval and must match."
    )


def test_plot_scripts_do_not_hardcode_the_txt_station_extension():
    """BP8 station files are .dat; a .txt-only glob silently matches nothing.

    scripts/plotDomainSweep.py looked for fltst_strk+000dp+000.txt, which no
    BP8 run produces, so it reported "no cases" for every BP8 sweep -- a result
    indistinguishable from having no runs at all. Nothing failed; the sweep just
    quietly produced an empty plot.
    """
    import glob as _glob
    import os
    import re as _re
    from conftest import ROOT, read

    bad = []
    for path in _glob.glob(os.path.join(str(ROOT), "scripts", "*.py")):
        src = read(os.path.join("scripts", os.path.basename(path)))
        if _re.search(r'fltst_strk[^"\']*\.txt', src):
            bad.append(os.path.basename(path))
    assert not bad, (
        f"{bad} look for BP8 station files with a .txt extension. BP8 writes "
        f".dat; match on 'fltst_strk...*' so both are accepted."
    )
