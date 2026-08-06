"""The BP8 initial condition must follow the benchmark, eq. (30).

    theta(x2, x3, 0) = D_RS / V_init

"The initial state variable is chosen at steady state with slip rate V_init
over the entire fault" -- SEAS BP8-QD-GS/PW, section 3.

This is prescribed, not derived. The guard exists because it is tempting to
"fix" it, and I did exactly that: eq. (30) together with Table 1 over-determines
the initial condition, and I resolved the contradiction by deriving theta from
tau_init instead. That silently made our submission non-conformant. Conformance
beats local repair -- a benchmark comparison is meaningless if each participant
patches the specification differently.

The over-determination, for the record:

    Table 1   tau_init = 14.6 MPa,  V_init = 1e-12 m/s,  sigma_bar_0 = 25 MPa
    eq. (30)  theta_0  = D_RS/V_init = 5.0e8 s
    eq. (13)  f(V_init, theta_0) = 0.51707  ->  tau = 12.93 MPa

12.93 != 14.6, so the fault does not start in equilibrium and the solver's
first step is V = 6.54e-11 m/s, not V_init. That discontinuity is a property of
the benchmark as written. It has been raised with the authors; if they revise
the specification, change these tests deliberately and say so in the commit.
"""

import math

import pytest

from conftest import load_case_params

BP8_CASES = ["bp8.qdc.gs.10", "test.bp8.qdc"]


def regularized_friction(a, b, Dc, v0, r0, V, theta):
    """BP8 eq. (13)."""
    psi = r0 + b * math.log(v0 * theta / Dc)
    return a * math.asinh(V / (2.0 * v0) * math.exp(psi / a))


@pytest.mark.parametrize("case", BP8_CASES)
def test_initial_state_follows_equation_30(case):
    """theta_0 must be exactly D_RS/V_init, everywhere on the fault."""
    par = load_case_params(case)
    expected = par.fric_rsf_Dc / par.init_slip_rate

    ofv = par.on_fault_vars
    nz, nx = ofv.shape[0], ofv.shape[1]
    probes = [(0, 0), (0, nx - 1), (nz - 1, 0), (nz - 1, nx - 1), (nz // 2, nx // 2)]

    for iz, ix in probes:
        theta = ofv[iz, ix, 20]
        assert theta == pytest.approx(expected, rel=1e-12), (
            f"{case}: theta_0 at ({iz},{ix}) is {theta:.6e} s, but BP8 eq. (30) "
            f"requires D_RS/V_init = {expected:.6e} s.\n"
            f"Do not derive theta_0 from tau_init. Eq. (30) prescribes it, and "
            f"deriving it instead makes the run non-conformant."
        )


@pytest.mark.parametrize("case", BP8_CASES)
def test_table_1_parameters_are_unmodified(case):
    """Table 1 values must reach the solver as written."""
    par = load_case_params(case)
    ofv = par.on_fault_vars

    assert par.fric_rsf_a == 0.016
    assert par.fric_rsf_b == 0.010
    assert par.fric_rsf_Dc == 0.5e-3
    assert par.fric_rsf_v0 == 1e-6
    assert par.fric_rsf_r0 == 0.6
    assert par.init_slip_rate == 1e-12
    assert abs(par.init_norm) == 25.0e6
    assert par.init_shear == 14.6e6

    assert ofv[0, 0, 8] == pytest.approx(14.6e6), "tau^0 must reach the solver as 14.6 MPa"
    assert abs(ofv[0, 0, 7]) == pytest.approx(25.0e6), "sigma_bar_0 must be 25 MPa"


@pytest.mark.parametrize("case", BP8_CASES)
def test_the_known_specification_inconsistency_is_unchanged(case):
    """Pin the size of the eq. (30) / Table 1 contradiction.

    Not a defect in our code -- a documented property of the benchmark. Pinned
    so that if it ever changes, someone looks at why. If the authors revise the
    specification this test should be updated, not deleted.
    """
    par = load_case_params(case)
    ofv = par.on_fault_vars
    theta0 = ofv[0, 0, 20]
    sigma_bar = abs(ofv[0, 0, 7])
    tau0 = ofv[0, 0, 8]

    f = regularized_friction(
        par.fric_rsf_a, par.fric_rsf_b, par.fric_rsf_Dc,
        par.fric_rsf_v0, par.fric_rsf_r0, par.init_slip_rate, theta0,
    )
    tau_from_state = sigma_bar * f
    assert tau_from_state / 1e6 == pytest.approx(12.9277, abs=1e-3)

    # And the slip rate the solver therefore actually starts from.
    psi = par.fric_rsf_r0 + par.fric_rsf_b * math.log(
        par.fric_rsf_v0 * theta0 / par.fric_rsf_Dc
    )
    v_start = (
        2.0 * par.fric_rsf_v0
        * math.sinh(tau0 / (sigma_bar * par.fric_rsf_a))
        / math.exp(psi / par.fric_rsf_a)
    )
    assert v_start == pytest.approx(6.542e-11, rel=1e-3), (
        f"{case}: the initial slip rate implied by eq. (30) + Table 1 changed "
        f"to {v_start:.4e} m/s. Expected 6.542e-11 m/s. If the benchmark was "
        f"revised, update this test deliberately."
    )
