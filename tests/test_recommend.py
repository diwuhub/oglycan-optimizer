"""Tests for constrained next-method recommendation."""

import os

from oglycan.core import load_params, load_site_catalog
from oglycan.recommend import recommend_next_method


def _repo_root():
    return os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _load_etanercept_inputs():
    root = _repo_root()
    catalog = load_site_catalog(os.path.join(root, "sites", "etanercept.json"))
    params = load_params(os.path.join(root, "examples", "default_acquisition_params.json"))
    return catalog, params


def test_recommend_next_method_never_worse_without_constraints():
    catalog, params = _load_etanercept_inputs()
    result = recommend_next_method(catalog, params)
    assert result["improvement"] >= 0.0
    assert result["composite_after"] >= result["composite_before"]


def test_recommend_next_method_respects_max_runtime_constraint():
    catalog, params = _load_etanercept_inputs()
    result = recommend_next_method(catalog, params, {"max_runtime_min": 70})
    gradient = result["recommended_params"]["lc_gradient"]["gradient_time_min"]
    assert gradient <= 50


def test_recommend_next_method_reports_impossible_runtime_constraint():
    catalog, params = _load_etanercept_inputs()
    result = recommend_next_method(catalog, params, {"max_runtime_min": 30})
    assert result["constraint_violations"]
    assert any("max_runtime_min=30" in item for item in result["constraint_violations"])


def test_recommend_next_method_noops_for_tims():
    catalog, params = _load_etanercept_inputs()
    result = recommend_next_method(catalog, params, {"instrument_class": "tims"})
    assert result["recommended_params"] == params
    assert result["changes"] == []
    assert result["improvement"] == 0.0
    assert any("tims" in item for item in result["constraint_violations"])
