"""Tests for glycan compositions, masses, and SVG rendering."""

from __future__ import annotations

import json
from math import isclose
from pathlib import Path

import pytest

from oglycan.glycan import (
    CORE_TYPE_COMPOSITION,
    glycoform_mass,
    render_glycoform_svg,
)


def test_glycoform_mass_core1_galnac():
    assert isclose(glycoform_mass("Core1_GalNAc"), 203.079372, abs_tol=1e-4)


def test_glycoform_mass_core1():
    assert isclose(glycoform_mass("Core1"), 365.132196, abs_tol=1e-4)


def test_glycoform_mass_core1_sia():
    assert isclose(glycoform_mass("Core1_Sia"), 656.227612, abs_tol=1e-4)


def test_glycoform_mass_core2():
    assert isclose(glycoform_mass("Core2"), 568.211568, abs_tol=1e-4)


def test_glycoform_mass_core2_fuc():
    assert isclose(glycoform_mass("Core2_Fuc"), 714.269477, abs_tol=1e-4)


def test_glycoform_mass_core3():
    assert isclose(glycoform_mass("Core3"), 406.158744, abs_tol=1e-4)


def test_all_catalog_core_types_have_compositions():
    sites_dir = Path(__file__).resolve().parents[1] / "sites"
    seen_core_types = set()
    for path in sites_dir.glob("*.json"):
        catalog = json.loads(path.read_text(encoding="utf-8"))
        for site in catalog.get("sites", []):
            seen_core_types.update(site.get("core_types", []))
    assert seen_core_types.issubset(CORE_TYPE_COMPOSITION)


def test_render_glycoform_svg_for_every_defined_core_type():
    for core_type in CORE_TYPE_COMPOSITION:
        svg = render_glycoform_svg(core_type)
        assert svg.startswith("<svg")
        assert svg.endswith("</svg>")


def test_unknown_core_type_raises_value_error():
    with pytest.raises(ValueError):
        glycoform_mass("Unknown_Core")
    with pytest.raises(ValueError):
        render_glycoform_svg("Unknown_Core")
