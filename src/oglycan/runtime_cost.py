"""Runtime cost helpers for gradient-time sweeps."""

from __future__ import annotations

import copy

from .core import evaluate

DEFAULT_GRADIENT_TIMES = [45, 60, 90, 120, 150, 180, 240]


def runtime_minutes(
    params: dict,
    n_injections: int = 1,
    overhead_min: float = 20.0,
) -> float:
    """Approximate total runtime minutes."""
    gradient_time = params.get("lc_gradient", {}).get("gradient_time_min", 120)
    return float(n_injections) * (float(gradient_time) + float(overhead_min))


def pareto_gradient_sweep(
    site_catalog: dict,
    params: dict,
    gradient_times_min: list = None,
    overhead_min: float = 20.0,
) -> list[dict]:
    """Evaluate composite and runtime across a fixed gradient-time sweep."""
    sweep = DEFAULT_GRADIENT_TIMES if gradient_times_min is None else list(gradient_times_min)
    points = []

    for gradient_time in sorted(sweep):
        trial_params = copy.deepcopy(params)
        trial_params["lc_gradient"]["gradient_time_min"] = gradient_time
        points.append({
            "gradient_time_min": gradient_time,
            "composite": evaluate(site_catalog, trial_params)["composite_score"],
            "runtime_min": runtime_minutes(trial_params, overhead_min=overhead_min),
            "on_pareto_front": False,
        })

    for i, point in enumerate(points):
        dominated = False
        for j, other in enumerate(points):
            if i == j:
                continue
            if other["composite"] > point["composite"] and other["runtime_min"] < point["runtime_min"]:
                dominated = True
                break
        point["on_pareto_front"] = not dominated

    return points
