"""HCD glycan fragment sub-model.

Scores oxonium ion and Y-ion generation for glycan identification.

Key trade-off: higher SA% improves HCD fragments but hurts ETD.
CE (collision energy) optimum is ~28-32 NCE; too low = no fragments,
too high = over-fragmentation losing diagnostic ions.

Literature basis:
  - Nilsson, Mol Cell Proteomics 2012: HCD oxonium ion diagnostics
  - Reiding et al., Anal Chem 2018: NCE optimization for O-glycopeptides
"""

from ._math import gauss as _gauss


def score(frag: dict, glycan_db: dict) -> float:
    """Score HCD fragment generation quality.

    Parameters:
        frag: fragmentation parameters (ethcd_sa_percent, collision_energy_nce)
        glycan_db: glycan database config (include_core1, include_core2, etc.)

    Returns:
        float in [0, 1]. Higher = better glycan fragment identification.
    """
    sa = frag.get("ethcd_sa_percent", 25)
    ce = frag.get("collision_energy_nce", 30)

    # CE optimum at 30 NCE
    ce_score = _gauss(ce, 30.0, 6.0)

    # SA% boosts HCD: optimum for HCD component is ~33-35%
    sa_hcd = _gauss(sa, 33.0, 14.0)

    # Glycan database completeness
    db = 0.0
    if glycan_db.get("include_core1"):
        db += 0.30
    if glycan_db.get("include_core2"):
        db += 0.25
    if glycan_db.get("include_sialylated"):
        db += 0.20
    if glycan_db.get("include_fucosylated"):
        db += 0.15
    if glycan_db.get("include_core3"):
        db += 0.05

    # Max glycan size: too small misses complex forms, too large adds FP
    max_sz = glycan_db.get("max_glycan_size", 4)
    sz_score = _gauss(max_sz, 5, 1.5)
    db = db * (0.5 + 0.5 * sz_score)

    return 0.35 * ce_score + 0.30 * sa_hcd + 0.35 * db
