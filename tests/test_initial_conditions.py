"""The BP8 initial condition must start the fault at a uniform V_init.

Section 3 of the benchmark description (2026-08-06 revision):

    "The initial state and reference shear traction on the fault is chosen so
     that the model can start with a uniform fault slip rate, given by
     V = [V_init, V_zero]"                                          -- eq. (27)
    "The initial state variable is chosen at steady state with slip rate V_init
     over the entire fault, namely theta(x2, x3, 0) = D_RS/V_init"  -- eq. (29)

Both `theta_0` and `tau^0` are *chosen* to make the fault start at `V_init`.
Eq. (28) only fixes the direction of the traction, `tau^0 = tau^0 * V/|V|`.
Eq. (29) pins `theta_0`, so the scalar `tau^0` is determined -- 12.9277 MPa.

Table 1 lists `tau_init` = 14.6 MPa, which contradicts that. The two cannot both
hold: 14.6 MPa with eq. (29)'s `theta_0` starts the fault at 6.542e-11 m/s, 65x
`V_init`, which is precisely what eq. (27) forbids. The prose is taken as
authoritative and Table 1's entry as stale -- the description is visibly derived
from BP6 and still carries BP6 text elsewhere ("(for BP6-A/S) state" in section
3, "for BP6-C set to zero" in section 4.1). Raised with the authors; if they
confirm 14.6 MPa is intended, change these tests deliberately and say so.

The derivation uses the same `shear_steady_state()` helper as the BP5 and BP7
compsets. Deviating from that convention is what produced the defect this guard
now prevents: `tau^0` was transcribed from Table 1 instead of derived, and every
BP8 result was ~1.8x too much slip.
"""

import math

import pytest

from conftest import load_case_params

BP8_CASES = ["bp8.qdc.gs.10", "test.bp8.qdc"]

# Slip rate spans decades; compare in log space. 1 % of a decade is far tighter
# than anything physically meaningful and far looser than round-off.
LOG_TOL = 0.01


def slip_rate_from(a, b, Dc, v0, r0, sigma_bar, tau, theta):
    """Invert eq. (12) for V given the state and stress the case file builds."""
    psi = r0 + b * math.log(v0 * theta / Dc)
    return 2.0 * v0 * math.sinh(tau / (sigma_bar * a)) / math.exp(psi / a)


@pytest.mark.parametrize("case", BP8_CASES)
def test_initial_state_follows_equation_29(case):
    """theta_0 must be exactly D_RS/V_init, everywhere on the fault."""
    par = load_case_params(case)
    expected = par.fric_rsf_Dc / par.init_slip_rate

    ofv = par.on_fault_vars
    nz, nx = ofv.shape[0], ofv.shape[1]
    probes = [(0, 0), (0, nx - 1), (nz - 1, 0), (nz - 1, nx - 1), (nz // 2, nx // 2)]

    for iz, ix in probes:
        theta = ofv[iz, ix, 20]
        assert theta == pytest.approx(expected, rel=1e-12), (
            f"{case}: theta_0 at ({iz},{ix}) is {theta:.6e} s, but eq. (29) "
            f"requires D_RS/V_init = {expected:.6e} s."
        )


@pytest.mark.parametrize("case", BP8_CASES)
def test_fault_starts_at_uniform_v_init(case):
    """The whole point of section 3: V(0) = V_init everywhere."""
    par = load_case_params(case)
    ofv = par.on_fault_vars
    nz, nx = ofv.shape[0], ofv.shape[1]
    probes = [(0, 0), (0, nx - 1), (nz - 1, 0), (nz - 1, nx - 1), (nz // 2, nx // 2)]

    for iz, ix in probes:
        v = slip_rate_from(
            par.fric_rsf_a, par.fric_rsf_b, par.fric_rsf_Dc,
            par.fric_rsf_v0, par.fric_rsf_r0,
            abs(ofv[iz, ix, 7]), ofv[iz, ix, 8], ofv[iz, ix, 20],
        )
        err = abs(math.log10(v) - math.log10(par.init_slip_rate))
        assert err < LOG_TOL, (
            f"{case}: the initial condition at ({iz},{ix}) starts the fault at "
            f"{v:.4e} m/s, not V_init = {par.init_slip_rate:.4e} m/s "
            f"({v/par.init_slip_rate:.1f}x off).\n"
            f"  tau^0   = {ofv[iz,ix,8]/1e6:.6f} MPa\n"
            f"  theta_0 = {ofv[iz,ix,20]:.4e} s\n"
            f"Section 3 requires a uniform start at V_init. Derive tau^0 with "
            f"shear_steady_state(), as the BP5 and BP7 compsets do -- do not "
            f"transcribe Table 1's tau_init."
        )


@pytest.mark.parametrize("case", BP8_CASES)
def test_derived_tau0_and_the_table_1_discrepancy(case):
    """Pin the derived tau^0 and the size of the contradiction with Table 1."""
    par = load_case_params(case)
    tau0 = par.on_fault_vars[0, 0, 8]
    assert tau0 / 1e6 == pytest.approx(12.9277, abs=1e-3), (
        f"{case}: derived tau^0 is {tau0/1e6:.6f} MPa, expected 12.9277 MPa."
    )
    # Table 1's value is retained but must not be the one that reaches the solver.
    assert hasattr(par, "init_shear_table1")
    assert par.init_shear_table1 == 14.6e6
    gap = (par.init_shear_table1 - tau0) / 1e6
    assert gap == pytest.approx(1.6723, abs=1e-3), (
        f"{case}: the Table 1 / section 3 gap moved to {gap:.4f} MPa "
        f"(was 1.6723). If the description was revised, update deliberately."
    )


@pytest.mark.parametrize("case", BP8_CASES)
def test_table_1_parameters_are_unmodified(case):
    par = load_case_params(case)
    ofv = par.on_fault_vars
    assert par.fric_rsf_a == 0.016
    assert par.fric_rsf_b == 0.010
    assert par.fric_rsf_Dc == 0.5e-3
    assert par.fric_rsf_v0 == 1e-6
    assert par.fric_rsf_r0 == 0.6
    assert par.init_slip_rate == 1e-12
    assert abs(par.init_norm) == 25.0e6
    assert abs(ofv[0, 0, 7]) == pytest.approx(25.0e6), "sigma_bar_0 must be 25 MPa"
