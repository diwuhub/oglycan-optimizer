"""Per-site localization confidence scoring.

Takes sub-model scores and applies per-site difficulty scaling and
glycoform diversity factors to produce individual site confidence values.
"""


def score_site(site: dict, etd: float, hcd: float, op_bonus: float,
               search: float, lc: float, ms: float) -> float:
    """Localization confidence for a single O-glycosite.

    Args:
        site: site dict with 'difficulty' (0-1) and 'core_types' (list).
        etd: ETD efficiency sub-model score.
        hcd: HCD fragment sub-model score.
        op_bonus: OpeRATOR enzyme bonus (0-0.35).
        search: search quality sub-model score.
        lc: LC separation sub-model score.
        ms: MS acquisition sub-model score.

    Returns:
        float in [0, 1]. Localization confidence for this site.
    """
    difficulty = site["difficulty"]
    n_forms = len(site["core_types"])

    # Fixed weights sum to 0.85; op_bonus (0-0.35) fills the headroom
    base = (0.25 * etd + 0.20 * hcd + 0.15 * ms + 0.15 * search
            + 0.10 * lc + op_bonus)

    # Why: harder sites must not score higher than easier sites under identical params.
    base *= (1.0 - 0.35 * difficulty)

    form_factor = 0.80 + 0.07 * min(n_forms, 3)
    base *= form_factor

    return max(0.0, min(1.0, base))
