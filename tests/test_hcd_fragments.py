"""Tests for HCD fragment sub-model."""
from oglycan.models.hcd_fragments import score

FULL_DB = {
    "include_core1": True, "include_core2": True,
    "include_sialylated": True, "include_fucosylated": True,
    "include_core3": True, "max_glycan_size": 5,
}


def test_ce_peak_at_30():
    frag_30 = {"collision_energy_nce": 30, "ethcd_sa_percent": 33}
    frag_20 = {**frag_30, "collision_energy_nce": 20}
    frag_45 = {**frag_30, "collision_energy_nce": 45}
    assert score(frag_30, FULL_DB) > score(frag_20, FULL_DB)
    assert score(frag_30, FULL_DB) > score(frag_45, FULL_DB)


def test_glycan_db_completeness():
    frag = {"collision_energy_nce": 30, "ethcd_sa_percent": 33}
    empty_db = {"include_core1": False, "include_core2": False, "max_glycan_size": 4}
    assert score(frag, FULL_DB) > score(frag, empty_db)


def test_score_in_range():
    for ce in [10, 20, 30, 40, 50]:
        assert 0.0 <= score({"collision_energy_nce": ce}, FULL_DB) <= 1.0
