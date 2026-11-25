"""Physical-behaviour invariants that must hold for any correct run of these
compset, independent of resolution or step count -- the class of check this
project's e2e tier was missing. All three multi-fault bugs recorded in
PROJECT_RULES.md (rule 12) produced a run that *completed* and printed
plausible-looking numbers; none violated a numeric regression lock on the
first run that hit them, because there was no prior lock to violate. A
physical invariant does not need a frozen prior value to catch a fault with
zero nodes, an all-zero slab, or effective normal stress that has gone
tensile.

These read the committed data/*/ files directly -- no build, no
solve, milliseconds -- so they run on every save and still catch a
regenerated gold that violates physics, not just a gold that moved.
"""

import csv
import json

import numpy as np
import pytest

from conftest import ROOT

def _field_references():
    """{label: fault-plane CSV} for every reference in the tree, discovered.

    A hardcoded list is the last residue of the per-benchmark harness: it went
    stale the moment BP5-dip90, BP7 and the step-over were retired. Ask the
    tree instead, so adding or removing a reference needs no edit here.
    """
    out = {}
    base = ROOT / "data"
    if not base.is_dir():
        return out
    for bench in sorted(p for p in base.iterdir() if p.is_dir()):
        for cand in [bench] + sorted(d for d in bench.iterdir()
                                     if d.is_dir()
                                     and d.name not in ("plots", "archive")):
            csvs = sorted(f for f in cand.glob("fault.*.csv")
                          if not f.name.endswith(".r.csv"))
            if csvs:
                out[str(cand.relative_to(base))] = csvs[-1]
    return out


FIELD_GOLD = _field_references()

BP8_SNAPSHOT = ROOT / "data" / "test.bp8.qdc.gs.10" /  "fault.05301.csv"
# There was never a CSV step-over gold: the two-fault case is BP1002, and its
# snapshots are netCDF, not CSV. The constant sat here unused and unreachable,
# so nothing noticed.
STEPOVER_SNAPSHOT = ROOT / "data" / "bp1002.qdc.2500" / "cycle0" / "fault.00001.nc"

BP8_GOLD_DIR = ROOT / "data" / "test.bp8.qdc.gs.10" 
BP8_STATIONS = [f"{s:+04d}dp{d:+04d}" for s in (-200, 0, 200) for d in (-200, 0, 200)]
T_OFF_S = 100 * 3600.0  # BP8 Table 1: injection turns off at 100 hours


def _read_csv(path):
    with open(path) as f:
        r = csv.DictReader(f)
        rows = list(r)
    return rows


def _col(rows, name):
    return np.array([float(r[name]) for r in rows])


@pytest.mark.parametrize("bench", list(FIELD_GOLD))
def test_effective_normal_stress_stays_compressive(bench):
    """Negative is compressive in this code's convention (par.init_norm is
    negative). A positive (tensile) effective normal stress means the fault
    has been pulled open, which none of these compset' physics produces."""
    rows = _read_csv(FIELD_GOLD[bench])
    en = _col(rows, "effective_normal")
    assert np.all(en <= 0.0), (
        f"{bench}: {int((en > 0).sum())} of {len(en)} fault nodes have "
        f"positive (tensile) effective normal stress, max = {en.max():.3e} Pa"
    )

def test_bp8_effective_normal_stress_stays_compressive():
    rows = _read_csv(BP8_SNAPSHOT)
    en = _col(rows, "effective_normal")
    assert np.all(en <= 0.0)


@pytest.mark.parametrize("bench", list(FIELD_GOLD))
def test_state_variable_stays_positive_and_finite(bench):
    """The rate-and-state state variable is a positive, finite quantity by
    construction (it has units of time in the aging law). A negative,
    zero, NaN or Inf state means the evolution law was fed a bad step or an
    uninitialised node reached the output -- exactly bug #2's signature."""
    rows = _read_csv(FIELD_GOLD[bench])
    st = _col(rows, "state_variable")
    assert np.all(np.isfinite(st)), f"{bench}: non-finite state_variable present"
    assert np.all(st > 0.0), f"{bench}: {int((st <= 0).sum())} nodes have state_variable <= 0"

def test_bp8_state_variable_stays_positive_and_finite():
    rows = _read_csv(BP8_SNAPSHOT)
    st = _col(rows, "state_variable")
    assert np.all(np.isfinite(st))
    assert np.all(st > 0.0)


# --- BP8 station time series: slip monotonicity and pore-pressure shape ----

def _station_rows(name):
    path = BP8_GOLD_DIR / f"fltst_strk{name}.csv"
    return _read_csv(path)


@pytest.mark.parametrize("station", BP8_STATIONS)
def test_bp8_slip_is_non_decreasing(station):
    """Cumulative slip cannot go backwards. A decrease would mean either a
    sign error in the slip accumulator or state overwritten by a neighbour
    (this project's accumulator-aliasing bug reached exactly this kind of
    quantity)."""
    rows = _station_rows(station)
    slip = _col(rows, "slip_2")
    diffs = np.diff(slip)
    bad = np.where(diffs < -1e-12)[0]  # float noise tolerance, not physical slack
    assert bad.size == 0, (
        f"{station}: slip_2 decreases at {bad.size} of {len(diffs)} steps, "
        f"first at row {bad[0]+1} ({slip[bad[0]]:.3e} -> {slip[bad[0]+1]:.3e} m)"
    )


def test_bp8_late_time_pressure_is_nearly_spatially_uniform():
    """Omega_f is closed (zero flux on its boundary; porepressure.f90's
    comment on pfFx/pfFz) and injection stops at t_off = 100 h, so every
    drop of injected fluid stays inside Omega_f and, given enough time,
    diffuses toward a spatially uniform pressure -- there is nowhere else
    for it to go and nothing left driving a gradient. This is the discrete
    fingerprint of that conservation: at t = 30 d the 9 stations' late-time
    pressures already agree to <1 % of their mean. A source or boundary bug
    that leaked mass out of Omega_f, or mis-normalized weights (the
    integrated-source identity src_i * A_i summing to q_inj/(beta*phi*Lfwid)
    documented in porepressure.f90), would show up as either a shrinking
    total or a persistent spatial gradient here.
    """
    summary = json.loads((BP8_GOLD_DIR / "summary.json").read_text())
    p_late = np.array([summary[s]["p_late_MPa"] for s in BP8_STATIONS])
    spread = (p_late.max() - p_late.min()) / p_late.mean()
    assert spread < 0.05, (
        f"late-time pore pressure spread across the 9 stations is "
        f"{spread*100:.2f} % of the mean ({p_late.min()}-{p_late.max()} MPa); "
        f"a closed, source-off system should be nearly uniform by t = 30 d"
    )


def test_bp8_pore_pressure_peaks_at_the_injection_point():
    """The Gaussian source is centred at (x2, x3) = (0, 0); every other
    station is farther from the injection point and must see a lower peak
    pressure. A wrong source location or a station/coordinate mixup would
    show up here as some off-centre station out-peaking the centre."""
    summary = json.loads((BP8_GOLD_DIR / "summary.json").read_text())
    center = summary["+000dp+000"]["p_peak_MPa"]
    for station in BP8_STATIONS:
        if station == "+000dp+000":
            continue
        other = summary[station]["p_peak_MPa"]
        assert other <= center, (
            f"{station}: peak pore pressure {other} MPa exceeds the "
            f"injection-point peak {center} MPa"
        )


def test_bp8_pore_pressure_decays_after_injection_shutoff():
    """fluid_toff = 100 h (BP8 Table 1): after shutoff the source term is
    zero everywhere, so pressure at the injection point can only diffuse
    outward and fall, never rise further. Checked on the raw run's own
    values, not an assumed shape."""
    rows = _station_rows("+000dp+000")
    t = _col(rows, "t")
    p = _col(rows, "pore_pressure")
    after = p[t > T_OFF_S]
    assert after.size > 2, "not enough post-shutoff samples in the gold CSV"
    # Late-time pressure must be below the pre-shutoff peak, and the last
    # third of the post-shutoff series must be non-increasing on average
    # (allow the coarse 500-row subsample some local noise; the trend must
    # not still be rising).
    peak_by_shutoff = p[t <= T_OFF_S].max()
    assert after[-1] < peak_by_shutoff, (
        f"late-time pressure {after[-1]:.4f} MPa is not below the "
        f"pre-shutoff peak {peak_by_shutoff:.4f} MPa"
    )
    tail = after[len(after) // 2:]
    assert tail[-1] <= tail[0] + 1e-6, (
        f"pore pressure is still rising in the back half of the post-shutoff "
        f"series ({tail[0]:.4f} -> {tail[-1]:.4f} MPa)"
    )


def test_multifault_reference_seeds_exactly_one_fault():
    """BP1002's step-1 snapshot must show the seed on fault 0 and only fault 0.

    This is the assertion the multi-fault reference exists to support, and it
    was missing: the README described the discriminating signal in prose while
    no test read it. Every multi-fault bug this project has hit -- a fault's
    slab written to the wrong index, both faults reading fault 1's data, a
    fault silently meshed with zero nodes -- shows up here as the seed landing
    on the wrong fault, on both, or on neither, and none of them needs a frozen
    prior value to detect.
    """
    if not STEPOVER_SNAPSHOT.exists():
        pytest.skip(f"{STEPOVER_SNAPSHOT} not present")
    nc = pytest.importorskip("netCDF4")
    with nc.Dataset(STEPOVER_SNAPSHOT) as d:
        assert "nid_fault" in d.dimensions, (
            f"{STEPOVER_SNAPSHOT.name} has no nid_fault dimension; it is not a "
            "multi-fault snapshot, so it cannot demonstrate per-fault routing")
        nfault = len(d.dimensions["nid_fault"])
        assert nfault == 2, f"expected 2 faults, found {nfault}"
        v = np.asarray(d["slip_rate"][:], float)

    seeded = [int((v[i] > 1e-3).sum()) for i in range(nfault)]
    assert seeded[0] > 0, (
        "fault 0 carries the nucleation patch but no node on it exceeds "
        f"1 mm/s at step 1 (counts per fault: {seeded}). The seed was not "
        "routed to fault 0.")
    assert seeded[1] == 0, (
        f"fault 1 must start at the creep rate everywhere -- whether it "
        f"ruptures is the question the case asks -- but {seeded[1]} of its "
        f"nodes already exceed 1 mm/s at step 1 (counts per fault: {seeded}). "
        "Fault 1 is reading fault 0's on-fault input.")
    # The two faults must not be byte-identical: that is what an aliased read
    # looks like once both happen to carry plausible values.
    assert not np.array_equal(v[0], v[1]), (
        "both faults have identical slip-rate slabs at step 1, which is what "
        "a fault-aliased read of on_fault_vars_input.nc produces")


# --- full-cycle metric functions (tests/e2e/cases.py, board D6) -------------
# Synthetic series with hand-computable answers. The rate is log-linear between
# samples, so a crossing of 1e-3 between 1e-4 and 1e-2 is at the midpoint.

def _cases():
    import sys
    sys.path.insert(0, str(ROOT / "tests" / "e2e"))
    import cases
    return cases


def _series(onset_shift=0.0):
    """Rate 1e-5, 1e-5, 1e-4, 1e-2, 1e-1, 1e-2, 1e-4, 1e-5, 1e-5 at t = 0..8 s
    (sample 2 moved by onset_shift); moment column 2.0 inside the event."""
    t = [0.0, 1.0, 2.0 + onset_shift, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
    v = [1e-5, 1e-5, 1e-4, 1e-2, 1e-1, 1e-2, 1e-4, 1e-5, 1e-5]
    return [[ti, vi, 2.0 if vi >= 1e-3 else 0.0, 0.0] for ti, vi in zip(t, v)]


def _tagged(events):
    for e in events:
        e.update(station="s", slip=1.0)
    return events


def test_find_events_log_linear_crossings_are_exact():
    ev = _cases().find_events(_series())
    assert [(len(ev)), ev[0]["onset"], ev[0]["end"]] == [1, pytest.approx(2.5, abs=1e-14),
                                                          pytest.approx(5.5, abs=1e-14)]


def test_find_events_moment_is_trapezoid_over_the_crossing_window():
    # window t = 2.5, 3, 4, 5, 5.5; column 3 interpolated to 1, then 2, 2, 2, then 1:
    # 0.75 + 2 + 2 + 0.75 = 5.5
    assert _cases().find_events(_series())[0]["moment"] == pytest.approx(5.5, abs=1e-14)


def test_find_events_peak_is_window_maximum():
    assert _cases().find_events(_series())[0]["peak"] == 1e-1


def test_find_events_two_separate_events_are_two():
    g = _series() + [[9.0 + k, v, 0.0, 0.0] for k, v in enumerate([1e-5, 1e-2, 1e-1, 1e-5])]
    assert len(_cases().find_events(g)) == 2


def test_find_events_a_lone_sample_exactly_at_threshold_is_an_event():
    # "at or above": a single sample at 1e-3 is one event (a strict > finds none).
    g = [[0.0, 1e-5, 0, 0], [1.0, 1e-3, 0, 0], [2.0, 1e-5, 0, 0]]
    assert len(_cases().find_events(g)) == 1


def test_find_events_window_touching_first_sample_has_nan_onset_and_finite_end():
    g = [[0.0, 1e-1, 0, 0], [1.0, 1e-1, 0, 0], [2.0, 1e-5, 0, 0]]
    ev = _cases().find_events(g)[0]
    assert (ev["onset"] != ev["onset"]) and ev["end"] == ev["end"]


def test_log_crossing_nonpositive_rate_raises():
    with pytest.raises(ValueError):
        _cases().log_crossing(0.0, 0.0, 1.0, 1e-2, 1e-3)


def test_station_slip_over_interpolates_hypot_difference(tmp_path):
    p = tmp_path / "fltst_x.txt"
    p.write_text("0 0 0\n1 3 4\n2 6 8\n")
    assert _cases().station_slip_over(str(p), 0.5, 1.5) == pytest.approx(5.0, abs=1e-14)


def test_compare_full_cycle_metrics_identical_passes():
    c = _cases()
    ev = _tagged(c.find_events(_series()))
    assert c.compare_full_cycle_metrics(ev, [dict(e) for e in ev])[0]


def test_compare_full_cycle_metrics_shifted_onset_is_caught():
    c = _cases()
    ok, msg = c.compare_full_cycle_metrics(_tagged(c.find_events(_series())),
                                           _tagged(c.find_events(_series(onset_shift=1e-6))))
    assert not ok and "onset_s" in msg


def test_compare_full_cycle_metrics_nan_peak_fails():
    c = _cases()
    ref = _tagged(c.find_events(_series()))
    assert not c.compare_full_cycle_metrics(ref, [dict(e, peak=float("nan")) for e in ref])[0]


def test_compare_full_cycle_metrics_event_count_mismatch_fails():
    c = _cases()
    ref = _tagged(c.find_events(_series()))
    assert not c.compare_full_cycle_metrics(ref, ref + ref)[0]
