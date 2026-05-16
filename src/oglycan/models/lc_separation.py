"""LC separation sub-model.

Scores gradient time, slope, flow rate, and column selection for
O-glycopeptide chromatographic resolution.

Key trade-off: longer gradients help but with diminishing returns
and peak broadening past ~150 min.

Literature basis:
  - Neue, J Chromatogr A 2005: gradient optimization theory
  - Reiding et al., Anal Chem 2018: glycopeptide gradient recommendations
"""

from ._math import gauss as _gauss


def score(lc: dict) -> float:
    """Score LC separation quality.

    Parameters (from lc dict):
        gradient_time_min: gradient length in minutes (default 120)
        gradient_start_pct_b: starting %B (default 2)
        gradient_end_pct_b: ending %B (default 40)
        flow_rate_uL_min: flow rate in uL/min (default 300)
        column: column descriptor string (default "")

    Returns:
        float in [0, 1]. Higher = better chromatographic separation.
    """
    grad_time = lc.get("gradient_time_min", 120)
    start_b = lc.get("gradient_start_pct_b", 2)
    end_b = lc.get("gradient_end_pct_b", 40)
    flow = lc.get("flow_rate_uL_min", 300)

    time_score = _gauss(grad_time, 120.0, 50.0)
    if grad_time < 45:
        time_score *= 0.5

    slope = (end_b - start_b) / max(1, grad_time)
    slope_score = _gauss(slope, 0.30, 0.08)

    start_score = _gauss(start_b, 2.0, 2.0)
    end_score = _gauss(end_b, 40.0, 8.0)
    flow_score = _gauss(flow, 300.0, 80.0)

    column = lc.get("column", "").lower()
    if "c18" in column and "uhplc" in column:
        col_bonus = 1.05
    elif "c18" in column:
        col_bonus = 1.0
    elif "hilic" in column:
        col_bonus = 0.95
    else:
        col_bonus = 0.90

    raw = (0.25 * time_score + 0.25 * slope_score + 0.15 * start_score
           + 0.15 * end_score + 0.20 * flow_score)
    return min(1.0, raw * col_bonus)
