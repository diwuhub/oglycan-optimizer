"""Tests for MS acquisition sub-model."""
from oglycan.models.ms_acquisition import score

def _default_ms(**overrides):
    base = {
        "resolution_ms1": 120000, "resolution_ms2": 30000,
        "agc_target_ms1": 1e6, "agc_target_ms2": 5e4,
        "max_injection_ms1": 50, "max_injection_ms2": 200,
        "dynamic_exclusion_s": 30,
    }
    base.update(overrides)
    return base


def test_ms2_30k_beats_60k():
    assert score(_default_ms(resolution_ms2=30000)) > score(_default_ms(resolution_ms2=60000))


def test_ms2_30k_beats_15k():
    assert score(_default_ms(resolution_ms2=30000)) > score(_default_ms(resolution_ms2=15000))


def test_dynamic_exclusion_30s_optimal():
    s_30 = score(_default_ms(dynamic_exclusion_s=30))
    assert s_30 > score(_default_ms(dynamic_exclusion_s=5))
    assert s_30 > score(_default_ms(dynamic_exclusion_s=120))


def test_score_in_range():
    for res2 in [15000, 30000, 60000, 120000]:
        assert 0.0 <= score(_default_ms(resolution_ms2=res2)) <= 1.0
