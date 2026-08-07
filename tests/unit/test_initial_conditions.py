"""The BP8 initial condition: Table 1's tau^0, and theta_0 derived from it.

The specification is over-determined. Eq. (12) ties tau, V and theta, so only
two may be prescribed, and section 3 prescribes three:

    Table 1   tau_init = 14.6 MPa
    eq. (27)  V(0) = V_init = 1e-12 m/s
    eq. (29)  theta_0 = D_RS/V_init = 5.0e8 s

Together, Table 1 and eq. (29) give V(0) = 6.542e-11 m/s, 65x what eq. (27)
asks for. Something has to go.

`taehoKim_ref` resolves it by keeping Table 1 and eq. (27) and letting eq. (29)
go, and both halves of that are directly visible on the CRESCENT DET comparison
(see reference/bp8/README.md): its shear_stress_2 starts at 14.6 MPa, and its
slip_rate_2 starts at log10 = -12. So theta_0 is derived: 4.0188e11 s, not
5.0e8 s.

That is what these cases now do, and it is also the closest of the three
readings to the reference -- 37.30 mm against ~38 mm at (0,0) and 19.56 against
~21 at (-200,0), with an edge/centre ratio of 0.524 against 0.553.

Do not "fix" theta_0 back to eq. (29): that reintroduces the 65x start. Do not
derive tau^0 from eq. (29) either -- that gives 12.9277 MPa, contradicts Table 1
and the reference, and was measurably the worst of the three (23.97 mm at the
centre). Both mistakes have been made in this project already.
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
def test_tau0_is_table_1(case):
    """tau^0 must be Table 1's 14.6 MPa, which is what the reference uses."""
    par = load_case_params(case)
    ofv = par.on_fault_vars
    nz, nx = ofv.shape[0], ofv.shape[1]
    for iz, ix in [(0, 0), (nz - 1, nx - 1), (nz // 2, nx // 2)]:
        assert ofv[iz, ix, 8] == pytest.approx(14.6e6), (
            f"{case}: tau^0 at ({iz},{ix}) is {ofv[iz,ix,8]/1e6:.4f} MPa, "
            f"expected Table 1's 14.6 MPa."
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
def test_derived_theta0_and_the_equation_29_discrepancy(case):
    """Pin the derived theta_0 and how far it sits from eq. (29)."""
    par = load_case_params(case)
    theta0 = par.on_fault_vars[0, 0, 20]
    assert theta0 == pytest.approx(4.0188e11, rel=1e-3), (
        f"{case}: derived theta_0 is {theta0:.4e} s, expected 4.0188e11 s."
    )
    eq29 = par.fric_rsf_Dc / par.init_slip_rate
    assert eq29 == pytest.approx(5.0e8, rel=1e-9)
    assert theta0 / eq29 == pytest.approx(803.8, rel=1e-2), (
        f"{case}: the eq. (29) discrepancy moved to {theta0/eq29:.1f}x "
        f"(was 803.8x). If the description was revised, update deliberately."
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
