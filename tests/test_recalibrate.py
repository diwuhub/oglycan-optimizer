"""Tests for difficulty recalibration from pilot evidence."""

import copy
import os

from oglycan.core import load_site_catalog
from oglycan.recalibrate import recalibrate_difficulty


def _repo_root():
    return os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _load_etanercept_catalog():
    root = _repo_root()
    return load_site_catalog(os.path.join(root, "sites", "etanercept.json"))


def _site_by_pos(catalog: dict, pos: int) -> dict:
    for site in catalog["sites"]:
        if site["pos"] == pos:
            return site
    raise AssertionError(f"missing site at position {pos}")


def test_recalibrate_difficulty_updates_bumps_and_adds_without_mutating_input():
    catalog = _load_etanercept_catalog()
    original = copy.deepcopy(catalog)
    pilot = {
        "source": "canonical",
        "glycoprotein_name": "Etanercept",
        "sites": [
            {
                "pos": 266,
                "aa": "T",
                "observed_localization": 0.3,
                "n_spectra": 4,
                "observed_glycoforms": ["Core1_Sia", "Core2"],
            },
            {
                "pos": 270,
                "aa": "S",
                "observed_localization": 0.85,
                "n_spectra": 2,
                "observed_glycoforms": ["Core1"],
            },
        ],
        "metadata": {"source_file": "synthetic.json"},
    }

    recalibrated = recalibrate_difficulty(catalog, pilot)

    t266 = _site_by_pos(recalibrated, 266)
    assert t266["difficulty"] == 0.8

    t237 = _site_by_pos(recalibrated, 237)
    assert t237["difficulty"] == 0.65
    assert any("did not observe" in note for note in t237["notes"])

    added = _site_by_pos(recalibrated, 270)
    assert added["provisional"] is True
    assert added["difficulty"] == 0.6
    assert added["core_types"] == ["Core1"]

    assert catalog == original
    assert _site_by_pos(catalog, 266)["difficulty"] == 0.9
    assert not any(site["pos"] == 270 for site in catalog["sites"])
    assert recalibrated["recalibration_metadata"]["n_sites_updated"] == 1
    assert recalibrated["recalibration_metadata"]["n_sites_unobserved"] == 12
    assert recalibrated["recalibration_metadata"]["n_sites_added"] == 1
