"""Tests for markdown method-card export."""

from datetime import datetime, timezone
import os
import re
import time

from oglycan.core import evaluate, load_params, load_site_catalog
from oglycan.method_card import build_method_card


def _repo_root():
    return os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def test_build_method_card_contains_expected_sections_and_values():
    root = _repo_root()
    catalog = load_site_catalog(os.path.join(root, "sites", "etanercept.json"))
    params = load_params(os.path.join(root, "examples", "default_acquisition_params.json"))
    evaluation = evaluate(catalog, params)

    start = time.perf_counter()
    card = build_method_card(catalog, params, evaluation)
    elapsed = time.perf_counter() - start

    assert elapsed < 1.0
    assert len(card) >= 500
    assert "Etanercept" in card
    assert f"{evaluation['composite_score']:.4f}" in card
    assert "OpeRATOR" in card
    assert re.search(r"\b0\.\d{3,4}\b", card)
    assert datetime.now(timezone.utc).strftime("%Y-%m-%d") in card
    assert "must be transcribed into the vendor method editor" in card
    assert "# " in card or "## " in card
