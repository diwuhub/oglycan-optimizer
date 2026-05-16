"""Tests for ETD efficiency sub-model."""
from oglycan.models.etd_efficiency import score


def test_peak_at_18_sa():
    frag_18 = {"ethcd_sa_percent": 18, "charge_state_min": 2, "charge_state_max": 6, "isolation_window_mz": 1.6}
    frag_10 = {**frag_18, "ethcd_sa_percent": 10}
    frag_30 = {**frag_18, "ethcd_sa_percent": 30}
    s18 = score(frag_18)
    assert s18 > score(frag_10)
    assert s18 > score(frag_30)


def test_charge_window_too_narrow():
    wide = {"ethcd_sa_percent": 18, "charge_state_min": 2, "charge_state_max": 6, "isolation_window_mz": 1.6}
    narrow = {**wide, "charge_state_min": 4, "charge_state_max": 5}
    assert score(wide) > score(narrow)


def test_isolation_window_peak():
    frag_opt = {"ethcd_sa_percent": 18, "charge_state_min": 2, "charge_state_max": 6, "isolation_window_mz": 1.6}
    frag_wide = {**frag_opt, "isolation_window_mz": 4.0}
    assert score(frag_opt) > score(frag_wide)


def test_score_in_range():
    for sa in [5, 18, 25, 40, 60]:
        s = score({"ethcd_sa_percent": sa})
        assert 0.0 <= s <= 1.0
