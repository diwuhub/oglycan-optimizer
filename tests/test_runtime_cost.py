"""Tests for runtime cost helpers and Pareto sweep."""

import os

from oglycan.core import load_params, load_site_catalog
from oglycan.runtime_cost import pareto_gradient_sweep, runtime_minutes


def _repo_root():
    return os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _load_etanercept_inputs():
    root = _repo_root()
    catalog = load_site_catalog(os.path.join(root, "sites", "etanercept.json"))
    params = load_params(os.path.join(root, "examples", "default_acquisition_params.json"))
    return catalog, params


def test_runtime_minutes_defaults_to_gradient_plus_overhead():
    _, params = _load_etanercept_inputs()
    assert runtime_minutes(params) == 140.0


def test_pareto_gradient_sweep_returns_sorted_points():
    catalog, params = _load_etanercept_inputs()
    points = pareto_gradient_sweep(catalog, params)
    assert len(points) == 7
    assert [point["gradient_time_min"] for point in points] == [45, 60, 90, 120, 150, 180, 240]


def test_pareto_gradient_sweep_marks_front_points():
    catalog, params = _load_etanercept_inputs()
    points = pareto_gradient_sweep(catalog, params)
    assert any(point["on_pareto_front"] for point in points)


def test_shorter_gradient_has_lower_runtime_than_longer_gradient():
    _, params = _load_etanercept_inputs()
    short_params = {
        **params,
        "lc_gradient": {**params["lc_gradient"], "gradient_time_min": 30},
    }
    long_params = {
        **params,
        "lc_gradient": {**params["lc_gradient"], "gradient_time_min": 240},
    }
    assert runtime_minutes(short_params) < runtime_minutes(long_params)
