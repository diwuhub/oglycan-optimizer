"""Tests for LC separation sub-model."""
from oglycan.models.lc_separation import score

def _default_lc(**overrides):
    base = {
        "gradient_time_min": 120, "gradient_start_pct_b": 2,
        "gradient_end_pct_b": 40, "flow_rate_uL_min": 300,
        "column": "C18_UHPLC",
    }
    base.update(overrides)
    return base


def test_120min_near_optimal():
    assert score(_default_lc(gradient_time_min=120)) > score(_default_lc(gradient_time_min=30))


def test_diminishing_returns_past_150():
    assert score(_default_lc(gradient_time_min=120)) >= score(_default_lc(gradient_time_min=200))


def test_c18_uhplc_bonus():
    assert score(_default_lc(column="C18_UHPLC")) > score(_default_lc(column=""))


def test_score_in_range():
    for t in [20, 60, 120, 180, 240]:
        assert 0.0 <= score(_default_lc(gradient_time_min=t)) <= 1.0
