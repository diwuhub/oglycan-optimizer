"""Benchmark tests for curated glycoprotein site catalogs."""

import copy
import os

import pytest

from oglycan.core import evaluate, load_params, load_site_catalog


def _repo_root():
    return os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _evaluate_catalog(catalog_name: str) -> dict:
    root = _repo_root()
    catalog = load_site_catalog(os.path.join(root, "sites", catalog_name))
    params = load_params(os.path.join(root, "examples", "default_acquisition_params.json"))
    return evaluate(catalog, params)


def _load_default_params() -> dict:
    return load_params(os.path.join(_repo_root(), "examples", "default_acquisition_params.json"))


def _set_nested(mapping: dict, path: tuple[str, ...], value) -> None:
    current = mapping
    for key in path[:-1]:
        current = current[key]
    current[path[-1]] = value


def _catalog_score(catalog_name: str, params: dict) -> float:
    root = _repo_root()
    catalog = load_site_catalog(os.path.join(root, "sites", catalog_name))
    return evaluate(catalog, params)["composite_score"]


CATALOG_GOLDENS = [
    ("epo.json", 1, 1, 0.9595),
    ("ctla4_ig.json", 4, 4, 0.9419),
    ("iga1_hinge.json", 6, 6, 0.9361),
    ("atacicept.json", 2, 2, 0.9396),
]
CATALOG_NAMES = [
    "etanercept.json",
    "epo.json",
    "ctla4_ig.json",
    "iga1_hinge.json",
    "atacicept.json",
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
NUDGE_PARAMS = {
    ("lc_gradient", "gradient_time_min"): 135,
    ("fragmentation", "ethcd_sa_percent"): 20,
    ("ms_acquisition", "resolution_ms2"): 60000,
}


@pytest.mark.parametrize(
    "catalog_name,total_sites,sites_localized,composite_score",
    CATALOG_GOLDENS,
)
def test_catalog_benchmarks(catalog_name, total_sites, sites_localized, composite_score):
    result = _evaluate_catalog(catalog_name)
    assert result["total_sites"] == total_sites
    assert result["sites_localized"] == sites_localized
    assert result["composite_score"] == pytest.approx(composite_score, abs=1e-4)


@pytest.mark.parametrize(
    "catalog_name",
    CATALOG_NAMES,
)
def test_catalogs_respond_to_stress_params(catalog_name):
    default_params = _load_default_params()
    stress_params = copy.deepcopy(default_params)
    for path, value in STRESS_PARAMS.items():
        _set_nested(stress_params, path, value)

    default_score = _catalog_score(catalog_name, default_params)
    stressed_score = _catalog_score(catalog_name, stress_params)

    assert stressed_score <= default_score - 0.15


@pytest.mark.parametrize("catalog_name", CATALOG_NAMES)
def test_catalogs_respond_to_nudge_params(catalog_name):
    default_params = _load_default_params()
    nudged_params = copy.deepcopy(default_params)
    for path, value in NUDGE_PARAMS.items():
        _set_nested(nudged_params, path, value)

    default_score = _catalog_score(catalog_name, default_params)
    nudged_score = _catalog_score(catalog_name, nudged_params)

    assert abs(nudged_score - default_score) <= 0.10
