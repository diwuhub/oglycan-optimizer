"""ETD backbone fragmentation sub-model.

Scores c/z ion generation efficiency for O-glycopeptide site localization.

Key trade-off: SA% (supplemental activation) helps HCD but suppresses ETD.
The ETD optimum is ~18% SA; above ~30% the radical-driven fragmentation
degrades sharply. EThcD compromise at ~25% already costs ETD performance.

Literature basis:
  - Riley & Coon, Anal Chem 2018: EThcD SA% optimization for glycopeptides
  - Reiding et al., Anal Chem 2018: charge-state dependence of ETD efficiency
"""

from ._math import gauss as _gauss


def score(frag: dict) -> float:
    """Score ETD efficiency from fragmentation parameters.

    Parameters (from frag dict):
        ethcd_sa_percent: supplemental activation percentage (default 25)
        charge_state_min: minimum precursor charge state (default 2)
        charge_state_max: maximum precursor charge state (default 6)
        isolation_window_mz: isolation window in Th (default 2.0)

    Returns:
        float in [0, 1]. Higher = better ETD backbone fragmentation.
    """
    sa = frag.get("ethcd_sa_percent", 25)
    z_min = frag.get("charge_state_min", 2)
    z_max = frag.get("charge_state_max", 6)
    iso = frag.get("isolation_window_mz", 2.0)

    # SA%: ETD optimum at 18%, sharp falloff above 30%
    sa_score = _gauss(sa, 18.0, 6.0)

    # Charge states: glycopeptides are typically z=3-5
    charge_score = 1.0
    if z_min > 2:
        charge_score -= 0.25 * (z_min - 2)
    if z_min > 3:
        charge_score -= 0.20  # missing z=3 is very bad
    if z_max < 5:
        charge_score -= 0.20 * (5 - z_max)
    if z_max > 8:
        charge_score -= 0.10 * (z_max - 8)
    if z_max - z_min < 2:
        charge_score -= 0.30  # too narrow a window
    charge_score = max(0.05, min(1.0, charge_score))

    # Isolation window: optimum ~1.6 Th for glycopeptides
    iso_score = _gauss(iso, 1.6, 0.6)

    return 0.45 * sa_score + 0.35 * charge_score + 0.20 * iso_score
