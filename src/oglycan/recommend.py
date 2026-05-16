"""Constrained next-method recommendation via single-pass coordinate descent."""

from __future__ import annotations

import copy

from .core import evaluate
from .sensitivity import get_nested, set_nested

SEARCH_SPACE = [
    ("fragmentation.ethcd_sa_percent", [12, 15, 18, 21, 25, 30]),
    ("fragmentation.collision_energy_nce", [24, 27, 30, 33, 36]),
    ("ms_acquisition.resolution_ms2", [15000, 30000, 60000, 120000]),
    ("ms_acquisition.dynamic_exclusion_s", [15, 30, 45, 60, 90]),
    ("lc_gradient.gradient_time_min", [45, 60, 90, 120, 150, 180, 240]),
    ("enzyme_preprocessing.digestion_time_hours", [1, 2, 4, 6, 8]),
]


def recommend_next_method(
    site_catalog: dict,
    current_params: dict,
    constraints: dict | None = None,
) -> dict:
    """
    Coordinate-wise grid search over a fixed set of high-impact parameters.
    """
    constraints = {} if constraints is None else dict(constraints)
    starting_params = copy.deepcopy(current_params)
    recommended_params = copy.deepcopy(current_params)
    composite_before = _composite(site_catalog, starting_params)
    constraint_violations = []

    instrument_class = constraints.get("instrument_class")
    if instrument_class == "tims":
        return {
            "recommended_params": recommended_params,
            "composite_before": composite_before,
            "composite_after": composite_before,
            "improvement": 0.0,
            "changes": [],
            "constraint_violations": [
                "instrument_class=tims is not supported in v0; returning current parameters unchanged."
            ],
        }
    if instrument_class not in (None, "orbitrap"):
        constraint_violations.append(
            f"instrument_class={instrument_class!r} is not recognized; assuming orbitrap behavior."
        )

    for path, values in SEARCH_SPACE:
        candidate_values = list(values)
        enforce_constraint = False

        if path == "lc_gradient.gradient_time_min" and constraints.get("max_runtime_min") is not None:
            max_runtime = float(constraints["max_runtime_min"])
            limit = max_runtime - 20.0
            candidate_values = [value for value in candidate_values if value <= limit]
            if not candidate_values:
                constraint_violations.append(
                    f"max_runtime_min={max_runtime:g} leaves no valid gradient_time_min after 20 min overhead."
                )
                continue
            enforce_constraint = True

        current_value = get_nested(recommended_params, path)
        best_value = current_value
        if enforce_constraint and current_value not in candidate_values:
            best_composite = float("-inf")
        else:
            best_composite = _composite(site_catalog, recommended_params)

        for value in candidate_values:
            trial_params = copy.deepcopy(recommended_params)
            set_nested(trial_params, path, value)
            trial_composite = _composite(site_catalog, trial_params)
            if trial_composite > best_composite:
                best_composite = trial_composite
                best_value = value

        set_nested(recommended_params, path, best_value)

    composite_after = _composite(site_catalog, recommended_params)
    min_composite = constraints.get("min_composite")
    if min_composite is not None and composite_after < float(min_composite):
        constraint_violations.append(
            f"min_composite={float(min_composite):.4f} is unreachable; best composite is {composite_after:.4f}."
        )

    return {
        "recommended_params": recommended_params,
        "composite_before": composite_before,
        "composite_after": composite_after,
        "improvement": round(composite_after - composite_before, 4),
        "changes": _collect_changes(starting_params, recommended_params),
        "constraint_violations": constraint_violations,
    }


def _collect_changes(before: dict, after: dict) -> list[dict]:
    changes = []
    for path, _ in SEARCH_SPACE:
        before_value = get_nested(before, path)
        after_value = get_nested(after, path)
        if before_value != after_value:
            changes.append(
                {"path": path, "before": before_value, "after": after_value}
            )
    return changes


def _composite(site_catalog: dict, params: dict) -> float:
    return float(evaluate(site_catalog, params)["composite_score"])
