"""Tests for signed local sensitivity analysis."""

import copy
import os

from oglycan.core import load_params, load_site_catalog
from oglycan.sensitivity import get_nested, sensitivity, set_nested


def _repo_root():
    return os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _load_etanercept_inputs():
    root = _repo_root()
    catalog = load_site_catalog(os.path.join(root, "sites", "etanercept.json"))
    params = load_params(os.path.join(root, "examples", "default_acquisition_params.json"))
    return catalog, params


def _entry_by_path(rows, path):
    for row in rows:
        if row["param_path"] == path:
            return row
    raise AssertionError(f"missing sensitivity row for {path}")


def test_sensitivity_returns_sorted_dict_rows():
    catalog, params = _load_etanercept_inputs()
    rows = sensitivity(catalog, params)
    assert rows
    assert all(isinstance(row, dict) for row in rows)
    assert rows == sorted(rows, key=lambda row: row["abs_effect"], reverse=True)


def test_use_sialexo_has_large_effect_with_explicit_sign_convention():
    catalog, params = _load_etanercept_inputs()
    default_rows = sensitivity(catalog, params)
    default_entry = _entry_by_path(default_rows, "enzyme_preprocessing.use_sialexo")
    assert default_entry["abs_effect"] >= 0.1
    assert default_entry["delta_up"] < 0.0

    disabled_params = copy.deepcopy(params)
    set_nested(disabled_params, "enzyme_preprocessing.use_sialexo", False)
    disabled_rows = sensitivity(catalog, disabled_params)
    disabled_entry = _entry_by_path(disabled_rows, "enzyme_preprocessing.use_sialexo")
    assert disabled_entry["delta_up"] > 0.1


def test_ethcd_sa_percent_is_top_five_numeric_mover():
    catalog, params = _load_etanercept_inputs()
    rows = sensitivity(catalog, params)
    numeric_rows = [row for row in rows if row["kind"] == "numeric"]
    top_paths = [row["param_path"] for row in numeric_rows[:5]]
    assert "fragmentation.ethcd_sa_percent" in top_paths


def test_flipping_boolean_twice_returns_original_value():
    _, params = _load_etanercept_inputs()
    path = "enzyme_preprocessing.use_pngasef"
    original = get_nested(params, path)
    trial = copy.deepcopy(params)
    set_nested(trial, path, not original)
    set_nested(trial, path, not get_nested(trial, path))
    assert get_nested(trial, path) == original
