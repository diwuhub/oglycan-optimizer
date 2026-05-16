"""Difficulty recalibration from pilot evidence."""

from __future__ import annotations

import copy
import json
from datetime import datetime

from .pilot import validate_pilot


def recalibrate_difficulty(
    site_catalog: dict,
    pilot_result: dict,
    learning_rate: float = 0.5,
) -> dict:
    """
    Re-estimate site difficulties from observed pilot localization outcomes.
    """
    errors = validate_pilot(pilot_result)
    if errors:
        raise ValueError("Invalid pilot result: " + "; ".join(errors))

    catalog_name = site_catalog.get("glycoprotein", {}).get("name", "")
    pilot_name = pilot_result.get("glycoprotein_name", "")
    if catalog_name.lower() != pilot_name.lower():
        raise ValueError(
            f"Pilot glycoprotein {pilot_name!r} does not match catalog {catalog_name!r}"
        )

    updated_catalog = copy.deepcopy(site_catalog)
    observed_sites = {}
    for site in pilot_result.get("sites", []):
        key = (site["pos"], site["aa"].upper())
        if site.get("n_spectra", 0) <= 0:
            continue
        observed_sites[key] = copy.deepcopy(site)

    n_sites_updated = 0
    n_sites_unobserved = 0

    for site in updated_catalog.get("sites", []):
        key = (site["pos"], site["aa"].upper())
        old_difficulty = float(site["difficulty"])
        pilot_site = observed_sites.pop(key, None)

        if pilot_site is None or pilot_site.get("observed_localization") is None:
            site["difficulty"] = _clamp_difficulty(old_difficulty + 0.05)
            _append_note(
                site,
                "Pilot run did not observe this site; difficulty increased by 0.05.",
            )
            n_sites_unobserved += 1
            continue

        observed_localization = float(pilot_site["observed_localization"])
        new_difficulty = ((1.0 - learning_rate) * old_difficulty) + (
            learning_rate * (1.0 - observed_localization)
        )
        site["difficulty"] = _clamp_difficulty(new_difficulty)
        n_sites_updated += 1

    n_sites_added = 0
    for pilot_site in sorted(observed_sites.values(), key=lambda item: (item["pos"], item["aa"])):
        observed_glycoforms = sorted(pilot_site.get("observed_glycoforms") or [])
        core_types = [form for form in observed_glycoforms if form != "Unknown"]
        notes = ["Added from pilot evidence outside the original catalog."]
        if "Unknown" in observed_glycoforms:
            notes.append(
                "Unknown pilot glycoforms were retained separately and excluded from core_types for scoring."
            )
        updated_catalog.setdefault("sites", []).append(
            {
                "pos": pilot_site["pos"],
                "aa": pilot_site["aa"],
                "core_types": core_types,
                "pilot_observed_glycoforms": observed_glycoforms,
                "difficulty": 0.6,
                "provisional": True,
                "notes": notes,
            }
        )
        n_sites_added += 1

    updated_catalog["recalibration_metadata"] = {
        "source": pilot_result["source"],
        "learning_rate": learning_rate,
        "n_sites_updated": n_sites_updated,
        "n_sites_unobserved": n_sites_unobserved,
        "n_sites_added": n_sites_added,
        "timestamp": datetime.now().isoformat(timespec="seconds"),
    }
    return updated_catalog


def write_catalog(catalog: dict, path: str) -> None:
    """Write a recalibrated site catalog as JSON."""
    with open(path, "w") as f:
        json.dump(catalog, f, indent=2)
        f.write("\n")


def _append_note(site: dict, note: str) -> None:
    notes = site.setdefault("notes", [])
    if note not in notes:
        notes.append(note)


def _clamp_difficulty(value: float) -> float:
    return max(0.10, min(0.98, round(float(value), 4)))
