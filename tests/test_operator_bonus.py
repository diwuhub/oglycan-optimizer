"""Tests for OpeRATOR enzyme preprocessing sub-model."""
from oglycan.models.operator_bonus import score


def test_sialexo_bonus():
    with_sialexo = score({"primary": "OpeRATOR", "use_sialexo": True, "digestion_time_hours": 4, "temperature_C": 37})
    without = score({"primary": "OpeRATOR", "use_sialexo": False, "digestion_time_hours": 4, "temperature_C": 37})
    assert with_sialexo - without > 0.10


def test_pngasef_additive():
    base = score({"primary": "OpeRATOR", "use_sialexo": True, "use_pngasef": False, "digestion_time_hours": 4, "temperature_C": 37})
    with_png = score({"primary": "OpeRATOR", "use_sialexo": True, "use_pngasef": True, "digestion_time_hours": 4, "temperature_C": 37})
    assert with_png > base


def test_no_enzyme_gives_zero():
    assert score({"primary": "trypsin", "use_sialexo": False}) < 0.01


def test_score_capped_at_035():
    s = score({"primary": "OpeRATOR", "use_sialexo": True, "use_pngasef": True, "digestion_time_hours": 4, "temperature_C": 37})
    assert s <= 0.35
