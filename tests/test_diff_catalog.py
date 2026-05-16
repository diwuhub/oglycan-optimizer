"""Unit tests for catalog diff computation."""

from oglycan.cli import _compute_catalog_diff


def _catalog(name: str, sites: list[dict], threshold: float = 0.75) -> dict:
    return {
        "glycoprotein": {"name": name},
        "sites": sites,
        "localization_threshold": threshold,
    }


def test_compute_catalog_diff_identical():
    catalog = _catalog(
        "Toy",
        [{"pos": 12, "aa": "T", "difficulty": 0.6, "core_types": ["Core1"]}],
    )

    diff = _compute_catalog_diff(catalog, catalog)

    assert set(diff) == {
        "catalog_name_old",
        "catalog_name_new",
        "sites_changed",
        "sites_added",
        "sites_removed",
        "threshold_old",
        "threshold_new",
    }
    assert diff["catalog_name_old"] == "Toy"
    assert diff["catalog_name_new"] == "Toy"
    assert diff["sites_changed"] == []
    assert diff["sites_added"] == []
    assert diff["sites_removed"] == []
    assert diff["threshold_old"] == 0.75
    assert diff["threshold_new"] == 0.75


def test_compute_catalog_diff_one_site_changed():
    old_catalog = _catalog(
        "Toy",
        [{"pos": 12, "aa": "T", "difficulty": 0.9, "core_types": ["Core1"]}],
    )
    new_catalog = _catalog(
        "Toy",
        [{"pos": 12, "aa": "T", "difficulty": 0.775, "core_types": ["Core1", "Core1_Sia"]}],
    )

    diff = _compute_catalog_diff(old_catalog, new_catalog)

    assert diff["sites_added"] == []
    assert diff["sites_removed"] == []
    assert diff["sites_changed"] == [
        {
            "pos": 12,
            "aa": "T",
            "difficulty_old": 0.9,
            "difficulty_new": 0.775,
            "core_types_old": ["Core1"],
            "core_types_new": ["Core1", "Core1_Sia"],
        }
    ]


def test_compute_catalog_diff_one_site_added_one_removed():
    old_catalog = _catalog(
        "Toy",
        [{"pos": 12, "aa": "T", "difficulty": 0.6, "core_types": ["Core1"]}],
    )
    new_catalog = _catalog(
        "Toy",
        [{"pos": 14, "aa": "S", "difficulty": 0.5, "core_types": ["Core1", "Core1_Sia"]}],
        threshold=0.8,
    )

    diff = _compute_catalog_diff(old_catalog, new_catalog)

    assert diff["sites_changed"] == []
    assert diff["sites_added"] == [
        {"pos": 14, "aa": "S", "difficulty": 0.5, "core_types": ["Core1", "Core1_Sia"]}
    ]
    assert diff["sites_removed"] == [
        {"pos": 12, "aa": "T", "difficulty": 0.6, "core_types": ["Core1"]}
    ]
    assert diff["threshold_old"] == 0.75
    assert diff["threshold_new"] == 0.8
