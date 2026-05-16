"""Import smoke test for the Streamlit app module."""

import importlib
import importlib.util

import pytest


def test_web_app_import_has_main_when_streamlit_installed():
    if importlib.util.find_spec("streamlit") is None:
        pytest.skip("streamlit not installed")
    mod = importlib.import_module("oglycan.web.app")
    assert hasattr(mod, "main")


def test_render_annotated_sequence_html_marks_catalog_sites():
    if importlib.util.find_spec("streamlit") is None:
        pytest.skip("streamlit not installed")

    mod = importlib.import_module("oglycan.web.app")
    catalog = {
        "sites": [
            {"pos": 2, "aa": "S", "core_types": ["Core1"], "difficulty": 0.5},
            {"pos": 6, "aa": "S", "core_types": ["Core1_Sia"], "difficulty": 0.7},
        ],
        "predicted_sites": [
            {"pos": 4, "aa": "T", "p_glycosite": 0.8, "source": "scan_st"},
        ],
        "localization_threshold": 0.75,
    }
    site_scores = {
        (2, "S"): {"confidence": 0.91},
        (6, "S"): {"confidence": 0.61},
    }

    rendered = mod.render_annotated_sequence_html("ASATASAT", catalog, site_scores)

    assert "oglycan-pass" in rendered
    assert "oglycan-fail" in rendered
    assert "oglycan-predicted" in rendered
    assert "oglycan-unmarked" in rendered
