"""Tests for bootstrap composite confidence intervals."""

import os

from oglycan.core import load_params, load_site_catalog
from oglycan.uncertainty import composite_with_ci


def _repo_root():
    return os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _load_etanercept_inputs():
    root = _repo_root()
    catalog = load_site_catalog(os.path.join(root, "sites", "etanercept.json"))
    params = load_params(os.path.join(root, "examples", "default_acquisition_params.json"))
    return catalog, params


def test_composite_with_ci_shape_and_bounds():
    catalog, params = _load_etanercept_inputs()
    result = composite_with_ci(catalog, params, n_samples=128, seed=7)
    assert set(result) == {"point", "ci_low", "ci_high", "n_samples", "samples"}
    assert result["n_samples"] == 128
    assert len(result["samples"]) == 128
    # Why: `point` is composite at unperturbed params (sits at every Gaussian
    # peak). Any perturbation can only lower composite, so samples form a
    # one-sided distribution below `point`; we require ci_high <= point.
    assert result["ci_low"] <= result["ci_high"] <= result["point"]


def test_composite_with_ci_reproducible_for_same_seed():
    catalog, params = _load_etanercept_inputs()
    first = composite_with_ci(catalog, params, n_samples=128, seed=11)
    second = composite_with_ci(catalog, params, n_samples=128, seed=11)
    assert first == second


def test_etanercept_ci_width_is_sane():
    catalog, params = _load_etanercept_inputs()
    result = composite_with_ci(catalog, params, seed=0)
    width = result["ci_high"] - result["ci_low"]
    # Why: empirical CI width for Etanercept defaults with n=1000 samples is
    # ~0.26; ceiling set at 0.40 gives margin for Monte Carlo noise without
    # masking a regression that balloons the bootstrap.
    assert 0.005 <= width <= 0.40
