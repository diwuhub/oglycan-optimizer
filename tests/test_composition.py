"""Composition test: full Etanercept evaluation with default params."""
import json
import os
from oglycan.core import evaluate, load_params, load_site_catalog

def _repo_root():
    return os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def test_default_etanercept_composite():
    root = _repo_root()
    catalog = load_site_catalog(os.path.join(root, "sites", "etanercept.json"))
    params = load_params(os.path.join(root, "examples", "default_acquisition_params.json"))
    result = evaluate(catalog, params)
    assert result["total_sites"] == 13
    assert result["composite_score"] >= 0.85
    assert result["sites_localized"] >= 8


def test_all_sub_model_scores_positive():
    root = _repo_root()
    catalog = load_site_catalog(os.path.join(root, "sites", "etanercept.json"))
    params = load_params(os.path.join(root, "examples", "default_acquisition_params.json"))
    result = evaluate(catalog, params)
    for name, val in result["sub_model_scores"].items():
        assert val > 0.0


def test_site_results_count():
    root = _repo_root()
    catalog = load_site_catalog(os.path.join(root, "sites", "etanercept.json"))
    params = load_params(os.path.join(root, "examples", "default_acquisition_params.json"))
    result = evaluate(catalog, params)
    assert len(result["site_results"]) == 13


def test_hcd_only_penalizes_etd():
    root = _repo_root()
    catalog = load_site_catalog(os.path.join(root, "sites", "etanercept.json"))
    params = load_params(os.path.join(root, "examples", "default_acquisition_params.json"))
    result_ethcd = evaluate(catalog, params)
    params_hcd = json.loads(json.dumps(params))
    params_hcd["fragmentation"]["mode"] = "HCD"
    result_hcd = evaluate(catalog, params_hcd)
    assert result_ethcd["composite_score"] > result_hcd["composite_score"]
