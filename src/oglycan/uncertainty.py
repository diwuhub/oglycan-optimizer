"""Bootstrap confidence intervals for composite score.

The v0 sub-models store Gaussian optima as hardcoded constants, so this module
does not patch model code at runtime. Instead, for each uncertain optimum we
perturb the corresponding parameter value by ``-delta`` where
``delta ~ N(0, mu_sigma)``, using the identity
``gauss(x - delta, mu, sigma) == gauss(x, mu + delta, sigma)``.

Limitation: both ETD and HCD score functions read the same
``fragmentation.ethcd_sa_percent`` field but use different literature optima.
We therefore sample both SA uncertainties independently and apply their implied
shifts to the same shared parameter value in a single evaluation, which means
the two SA optima cannot be perturbed independently inside one model call.
"""

from __future__ import annotations

import copy
import math
import random

from .core import evaluate

OPTIMA_UNCERTAINTY = {
    # fragmentation
    "ethcd_sa_percent_etd_mu": (18.0, 2.0),
    "ethcd_sa_percent_hcd_mu": (33.0, 4.7),
    "collision_energy_nce_mu": (30.0, 2.0),
    "isolation_window_mz_mu": (1.6, 0.2),
    # ms_acquisition
    "resolution_ms2_mu": (30000.0, 2700.0),
    "cycle_time_mu": (2.5, 0.33),
    "dynamic_exclusion_s_mu": (30.0, 4.0),
    # lc_separation
    "gradient_time_min_mu": (120.0, 17.0),
    "gradient_slope_mu": (0.30, 0.027),
    "flow_rate_mu": (300.0, 27.0),
    # operator_bonus
    "digestion_time_hours_mu": (4.0, 0.67),
    "temperature_C_mu": (37.0, 1.33),
    # search_quality
    "mass_tolerance_ms1_ppm_mu": (8.0, 0.83),
    "mass_tolerance_ms2_ppm_mu": (18.0, 1.67),
    "min_peptide_length_mu": (6.5, 0.5),
    "max_missed_cleavages_mu": (2.0, 0.27),
    "fdr_log10_mu": (-2.0, 0.13),
}

# Maps uncertain optima to params dict paths. ``None`` means the optimum is
# structurally implied by multiple parameters rather than directly settable.
OPTIMA_PARAM_PATHS = {
    "ethcd_sa_percent_etd_mu": ("fragmentation.ethcd_sa_percent", "linear"),
    "ethcd_sa_percent_hcd_mu": ("fragmentation.ethcd_sa_percent", "linear"),
    "collision_energy_nce_mu": ("fragmentation.collision_energy_nce", "linear"),
    "isolation_window_mz_mu": ("fragmentation.isolation_window_mz", "linear"),
    "resolution_ms2_mu": ("ms_acquisition.resolution_ms2", "linear"),
    "cycle_time_mu": (None, "derived"),
    "dynamic_exclusion_s_mu": ("ms_acquisition.dynamic_exclusion_s", "linear"),
    "gradient_time_min_mu": ("lc_gradient.gradient_time_min", "linear"),
    "gradient_slope_mu": (None, "derived"),
    "flow_rate_mu": ("lc_gradient.flow_rate_uL_min", "linear"),
    "digestion_time_hours_mu": ("enzyme_preprocessing.digestion_time_hours", "linear"),
    "temperature_C_mu": ("enzyme_preprocessing.temperature_C", "linear"),
    "mass_tolerance_ms1_ppm_mu": ("search_params.mass_tolerance_ms1_ppm", "linear"),
    "mass_tolerance_ms2_ppm_mu": ("search_params.mass_tolerance_ms2_ppm", "linear"),
    "min_peptide_length_mu": ("search_params.min_peptide_length", "linear"),
    "max_missed_cleavages_mu": ("search_params.max_missed_cleavages", "linear"),
    "fdr_log10_mu": ("search_params.fdr_threshold", "log10"),
}


def _get_nested(data: dict, path: str):
    value = data
    for part in path.split("."):
        value = value[part]
    return value


def _set_nested(data: dict, path: str, value) -> None:
    target = data
    parts = path.split(".")
    for part in parts[:-1]:
        target = target[part]
    target[parts[-1]] = value


def _percentile(sorted_values: list[float], pct: float) -> float:
    if not sorted_values:
        raise ValueError("percentile requires at least one value")
    if len(sorted_values) == 1:
        return sorted_values[0]
    index = (len(sorted_values) - 1) * pct
    lower = int(math.floor(index))
    upper = int(math.ceil(index))
    if lower == upper:
        return sorted_values[lower]
    weight = index - lower
    return sorted_values[lower] * (1.0 - weight) + sorted_values[upper] * weight


def composite_with_ci(
    site_catalog: dict,
    params: dict,
    n_samples: int = 1000,
    seed: int = 0,
) -> dict:
    """Return bootstrap composite score CI under uncertain Gaussian optima."""
    if n_samples <= 0:
        raise ValueError("n_samples must be > 0")

    rng = random.Random(seed)
    point = evaluate(site_catalog, params)["composite_score"]
    samples = []

    for _ in range(n_samples):
        sampled_params = copy.deepcopy(params)
        linear_shifts = {}

        for key, (_, mu_sigma) in OPTIMA_UNCERTAINTY.items():
            path, space = OPTIMA_PARAM_PATHS[key]
            if path is None:
                # cycle_time and gradient_slope are derived from multiple knobs.
                continue

            delta = rng.gauss(0.0, mu_sigma)
            if space == "linear":
                linear_shifts[path] = linear_shifts.get(path, 0.0) + delta
                continue

            if space == "log10":
                current = float(_get_nested(params, path))
                current_log10 = math.log10(max(1e-12, current))
                _set_nested(sampled_params, path, 10 ** (current_log10 - delta))
                continue

            raise ValueError(f"unsupported perturbation space for {key}: {space}")

        for path, total_delta in linear_shifts.items():
            current = _get_nested(params, path)
            _set_nested(sampled_params, path, current - total_delta)

        samples.append(evaluate(site_catalog, sampled_params)["composite_score"])

    ordered = sorted(samples)
    return {
        "point": point,
        "ci_low": round(_percentile(ordered, 0.025), 4),
        "ci_high": round(_percentile(ordered, 0.975), 4),
        "n_samples": n_samples,
        "samples": samples,
    }
