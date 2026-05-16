"""Signed local sensitivity analysis for composite score."""

from __future__ import annotations

import copy

from .core import evaluate

NUMERIC_PATHS = [
    "fragmentation.ethcd_sa_percent",
    "fragmentation.collision_energy_nce",
    "fragmentation.isolation_window_mz",
    "ms_acquisition.resolution_ms1",
    "ms_acquisition.resolution_ms2",
    "ms_acquisition.agc_target_ms1",
    "ms_acquisition.agc_target_ms2",
    "ms_acquisition.max_injection_ms1",
    "ms_acquisition.max_injection_ms2",
    "ms_acquisition.dynamic_exclusion_s",
    "lc_gradient.gradient_time_min",
    "lc_gradient.gradient_start_pct_b",
    "lc_gradient.gradient_end_pct_b",
    "lc_gradient.flow_rate_uL_min",
    "enzyme_preprocessing.digestion_time_hours",
    "enzyme_preprocessing.temperature_C",
    "search_params.mass_tolerance_ms1_ppm",
    "search_params.mass_tolerance_ms2_ppm",
    "search_params.min_peptide_length",
    "search_params.max_missed_cleavages",
    "search_params.fdr_threshold",
    "glycan_database.max_glycan_size",
]

BOOLEAN_PATHS = [
    "enzyme_preprocessing.use_sialexo",
    "enzyme_preprocessing.use_pngasef",
    "glycan_database.include_core1",
    "glycan_database.include_core2",
    "glycan_database.include_core3",
    "glycan_database.include_sialylated",
    "glycan_database.include_fucosylated",
]

# TODO: support discrete string perturbations for these fields in a later phase.
DISCRETE_PATHS = [
    "enzyme_preprocessing.primary",
    "fragmentation.mode",
    "lc_gradient.column",
]


def get_nested(data: dict, path: str):
    """Read a dotted path from a nested dict."""
    value = data
    for part in path.split("."):
        value = value[part]
    return value


def set_nested(data: dict, path: str, value) -> None:
    """Write a dotted path inside a nested dict."""
    target = data
    parts = path.split(".")
    for part in parts[:-1]:
        target = target[part]
    target[parts[-1]] = value


def _composite(site_catalog: dict, params: dict) -> float:
    return evaluate(site_catalog, params)["composite_score"]


def sensitivity(
    site_catalog: dict,
    params: dict,
    delta_pct: float = 10.0,
) -> list[dict]:
    """
    For each perturbable parameter, compute signed sensitivity near the current point.
    """
    current_composite = _composite(site_catalog, params)
    pct_scale = delta_pct / 100.0
    rows = []

    for path in NUMERIC_PATHS:
        current_value = get_nested(params, path)
        if current_value == 0:
            continue

        up_params = copy.deepcopy(params)
        down_params = copy.deepcopy(params)
        set_nested(up_params, path, current_value * (1.0 + pct_scale))
        set_nested(down_params, path, current_value * (1.0 - pct_scale))

        delta_up = round(_composite(site_catalog, up_params) - current_composite, 4)
        delta_down = round(_composite(site_catalog, down_params) - current_composite, 4)
        rows.append({
            "param_path": path,
            "current_value": current_value,
            "kind": "numeric",
            "delta_up": delta_up,
            "delta_down": delta_down,
            "abs_effect": round(abs(delta_up) + abs(delta_down), 4),
        })

    for path in BOOLEAN_PATHS:
        current_value = bool(get_nested(params, path))
        flipped_params = copy.deepcopy(params)
        set_nested(flipped_params, path, not current_value)
        delta_up = round(_composite(site_catalog, flipped_params) - current_composite, 4)
        rows.append({
            "param_path": path,
            "current_value": current_value,
            "kind": "boolean",
            "delta_up": delta_up,
            "delta_down": 0.0,
            "abs_effect": round(abs(delta_up), 4),
        })

    rows.sort(key=lambda item: item["abs_effect"], reverse=True)
    return rows
