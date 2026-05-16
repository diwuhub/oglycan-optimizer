"""Tests for search quality sub-model."""
from oglycan.models.search_quality import score

def _default_search(**overrides):
    base = {
        "mass_tolerance_ms1_ppm": 10, "mass_tolerance_ms2_ppm": 20,
        "min_peptide_length": 6, "max_missed_cleavages": 2, "fdr_threshold": 0.01,
    }
    base.update(overrides)
    return base


def test_ms1_tolerance_peak():
    s_8 = score(_default_search(mass_tolerance_ms1_ppm=8))
    assert s_8 > score(_default_search(mass_tolerance_ms1_ppm=3))
    assert s_8 > score(_default_search(mass_tolerance_ms1_ppm=25))


def test_ms2_tolerance_peak():
    s_18 = score(_default_search(mass_tolerance_ms2_ppm=18))
    assert s_18 > score(_default_search(mass_tolerance_ms2_ppm=5))
    assert s_18 > score(_default_search(mass_tolerance_ms2_ppm=40))


def test_fdr_001_optimal():
    s_01 = score(_default_search(fdr_threshold=0.01))
    assert s_01 > score(_default_search(fdr_threshold=0.001))
    assert s_01 > score(_default_search(fdr_threshold=0.1))


def test_score_in_range():
    for tol in [3, 8, 15, 25]:
        assert 0.0 <= score(_default_search(mass_tolerance_ms1_ppm=tol)) <= 1.0
