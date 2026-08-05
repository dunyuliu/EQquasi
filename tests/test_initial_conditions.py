"""The prescribed initial condition must be self-consistent.

tau^0, V_init and theta^0 are linked by the regularized rate-and-state law:

    tau = sigma_bar * a * asinh[ V/(2 V*) * exp(Psi/a) ],
    Psi = f* + b * ln(V* theta / Dc)

Any two of them fix the third. A case file that prescribes all three
independently is over-determined, and the solver will silently start from
whatever V the (tau^0, theta^0) pair implies -- not the V_init that is
declared, documented and reported.

This is not hypothetical. Through v1.4.5 the BP8 cases carried

    theta^0 = Dc/V_init      (steady state, copied from the BP5 case)
    tau^0   = 14.6 MPa       (Table 1)

which are incompatible: steady state at V_init = 1e-12 m/s is 12.93 MPa, so
the fault began 1.67 MPa weaker than specified and ran at 6.5e-11 m/s -- 65x
V_init -- from the first step. Every BP8 result before that fix slipped about
50 % too far.

The guard walks each BP8-family case, evaluates the friction law at the
initial condition the file actually builds, and requires the resulting slip
rate to match the declared V_init.
"""

import math
import os
import sys

import pytest

from conftest import ROOT, load_case_params

BP8_CASES = ["bp8.qdc.gs.10", "test.bp8.qdc"]

# Slip rate spans decades, so compare in log space. 1 % of a decade is far
# tighter than any physically meaningful drift and far looser than the
# round-off of a correct closed-form inversion.
LOG_TOL = 0.01


def slip_rate_from(a, b, Dc, v0, r0, sigma_bar, tau, theta):
    """Invert eq. (9) for V given the state and stress the case file sets."""
    psi = r0 + b * math.log(v0 * theta / Dc)
    return 2.0 * v0 * math.sinh(tau / (sigma_bar * a)) / math.exp(psi / a)


@pytest.mark.parametrize("case", BP8_CASES)
def test_initial_slip_rate_matches_declared_v_init(case):
    par = load_case_params(case)

    a = par.fric_rsf_a
    b = par.fric_rsf_b
    Dc = par.fric_rsf_Dc
    v0 = par.fric_rsf_v0
    r0 = par.fric_rsf_r0
    v_init = par.init_slip_rate

    ofv = par.on_fault_vars
    nz, nx = ofv.shape[0], ofv.shape[1]

    # Corners and centre: enough to catch a uniform error and a bad ramp.
    probes = [(0, 0), (0, nx - 1), (nz - 1, 0), (nz - 1, nx - 1), (nz // 2, nx // 2)]

    for iz, ix in probes:
        theta = ofv[iz, ix, 20]
        sigma_bar = abs(ofv[iz, ix, 7])
        tau = ofv[iz, ix, 8]

        assert theta > 0.0, f"{case}: non-positive theta^0 at ({iz},{ix})"

        v = slip_rate_from(a, b, Dc, v0, r0, sigma_bar, tau, theta)
        err = abs(math.log10(v) - math.log10(v_init))

        assert err < LOG_TOL, (
            f"{case}: initial condition is over-determined at ({iz},{ix}).\n"
            f"  tau^0    = {tau/1e6:.4f} MPa\n"
            f"  sigma^0  = {sigma_bar/1e6:.4f} MPa\n"
            f"  theta^0  = {theta:.4e} s\n"
            f"  implied V = {v:.4e} m/s\n"
            f"  declared V_init = {v_init:.4e} m/s  ({v/v_init:.1f}x off)\n"
            f"Prescribe two of (tau^0, V_init, theta^0) and derive the third."
        )


@pytest.mark.parametrize("case", BP8_CASES)
def test_state_is_derived_not_hardcoded_to_steady_state(case):
    """theta = Dc/V_init is only right when tau^0 is the steady-state stress.

    Catches the specific regression directly, so the failure message names the
    mistake even if someone loosens the tolerance above.
    """
    par = load_case_params(case)
    theta = par.on_fault_vars[0, 0, 20]
    steady = par.fric_rsf_Dc / par.init_slip_rate

    tau_steady = (
        abs(par.on_fault_vars[0, 0, 7])
        * par.fric_rsf_a
        * math.asinh(
            par.init_slip_rate
            / (2.0 * par.fric_rsf_v0)
            * math.exp(
                (par.fric_rsf_r0 + par.fric_rsf_b * math.log(par.fric_rsf_v0 / par.init_slip_rate))
                / par.fric_rsf_a
            )
        )
    )
    tau = par.on_fault_vars[0, 0, 8]

    # If tau^0 really is the steady-state stress, theta == Dc/V_init is correct
    # and there is nothing to flag.
    if abs(tau - tau_steady) / tau_steady < 1e-6:
        return

    assert abs(theta - steady) / steady > 1e-6, (
        f"{case}: theta^0 is hard-coded to Dc/V_init = {steady:.4e} s, but "
        f"tau^0 = {tau/1e6:.4f} MPa is not the steady-state stress "
        f"({tau_steady/1e6:.4f} MPa). Derive theta^0 from tau^0 and V_init."
    )
