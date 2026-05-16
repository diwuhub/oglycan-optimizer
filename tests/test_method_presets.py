"""Face-validity tests for named method presets."""

import copy
import json
import os

import pytest

from oglycan.core import evaluate, load_params, load_site_catalog


CATALOG_NAMES = [
    "etanercept.json",
    "epo.json",
    "ctla4_ig.json",
    "iga1_hinge.json",
    "atacicept.json",
]
PRESET_NAMES = [
    "fast_throughput",
    "high_resolution",
    "minimal_enzyme",
]
STRESS_PARAMS = {
    ("fragmentation", "mode"): "HCD",
    ("fragmentation", "ethcd_sa_percent"): 40,
    ("ms_acquisition", "resolution_ms2"): 15000,
    ("lc_gradient", "gradient_time_min"): 30,
    ("enzyme_preprocessing", "use_sialexo"): False,
    ("enzyme_preprocessing", "use_pngasef"): False,
    ("search_params", "fdr_threshold"): 0.1,
}


def _repo_root():
    return os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _set_nested(mapping: dict, path: tuple[str, ...], value) -> None:
    current = mapping
    for key in path[:-1]:
        current = current[key]
    current[path[-1]] = value


def _load_default_params() -> dict:
    return load_params(os.path.join(_repo_root(), "examples", "default_acquisition_params.json"))


def _load_preset_params(preset_name: str) -> dict:
    return load_params(
        os.path.join(_repo_root(), "examples", "method_presets", f"{preset_name}.json")
    )


def _load_stress_params() -> dict:
    stress_params = copy.deepcopy(_load_default_params())
    for path, value in STRESS_PARAMS.items():
        _set_nested(stress_params, path, value)
    return stress_params


def _catalog_score(catalog_name: str, params: dict) -> float:
    root = _repo_root()
    catalog = load_site_catalog(os.path.join(root, "sites", catalog_name))
    return evaluate(catalog, params)["composite_score"]


def test_preset_ranking_on_etanercept():
    catalog_name = "etanercept.json"
    default_score = _catalog_score(catalog_name, _load_default_params())
    stress_score = _catalog_score(catalog_name, _load_stress_params())
    preset_scores = {
        preset_name: _catalog_score(catalog_name, _load_preset_params(preset_name))
        for preset_name in PRESET_NAMES
    }

    assert default_score >= preset_scores["fast_throughput"]
    assert default_score >= preset_scores["high_resolution"]
    assert default_score >= preset_scores["minimal_enzyme"]
    assert preset_scores["fast_throughput"] >= stress_score
    assert preset_scores["high_resolution"] >= stress_score
    assert preset_scores["minimal_enzyme"] >= stress_score


@pytest.mark.parametrize("catalog_name", CATALOG_NAMES)
def test_presets_rank_below_default(catalog_name):
    default_score = _catalog_score(catalog_name, _load_default_params())
    for preset_name in PRESET_NAMES:
        preset_score = _catalog_score(catalog_name, _load_preset_params(preset_name))
        assert preset_score <= default_score


@pytest.mark.parametrize("catalog_name", CATALOG_NAMES)
def test_presets_rank_above_stress(catalog_name):
    stress_score = _catalog_score(catalog_name, _load_stress_params())
    for preset_name in PRESET_NAMES:
        preset_score = _catalog_score(catalog_name, _load_preset_params(preset_name))
        assert preset_score > stress_score


def test_preset_metadata_fields_present():
    root = _repo_root()
    for preset_name in PRESET_NAMES:
        with open(
            os.path.join(root, "examples", "method_presets", f"{preset_name}.json"),
            encoding="utf-8",
        ) as handle:
            preset = json.load(handle)

        metadata = preset["_preset_metadata"]
        for key in ("name", "scenario", "derived_from"):
            assert isinstance(metadata[key], str)
            assert metadata[key]
