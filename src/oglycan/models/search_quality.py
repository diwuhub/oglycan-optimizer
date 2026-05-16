"""Database search parameter sub-model.

Scores mass tolerance, FDR threshold, peptide length, and missed
cleavage settings for glycopeptide identification.

Literature basis:
  - Eng et al., J Proteome Res 2008: mass tolerance optimization
  - Riley et al., J Proteome Res 2020: EThcD-specific search parameters
"""

import math

from ._math import gauss as _gauss


def score(search: dict) -> float:
    """Score database search parameter quality.

    Parameters (from search dict):
        mass_tolerance_ms1_ppm: MS1 mass tolerance in ppm (default 10)
        mass_tolerance_ms2_ppm: MS2 mass tolerance in ppm (default 20)
        min_peptide_length: minimum peptide length (default 6)
        max_missed_cleavages: maximum missed cleavages (default 2)
        fdr_threshold: FDR threshold (default 0.01)

    Returns:
        float in [0, 1]. Higher = better search parameter configuration.
    """
    tol1 = search.get("mass_tolerance_ms1_ppm", 10)
    tol2 = search.get("mass_tolerance_ms2_ppm", 20)
    min_len = search.get("min_peptide_length", 6)
    missed = search.get("max_missed_cleavages", 2)
    fdr = search.get("fdr_threshold", 0.01)

    tol1_score = _gauss(tol1, 8.0, 2.5)
    tol2_score = _gauss(tol2, 18.0, 5.0)
    len_score = _gauss(min_len, 6.5, 1.5)
    mc_score = _gauss(missed, 2.0, 0.8)
    fdr_score = _gauss(math.log10(max(1e-6, fdr)), -2.0, 0.4)

    return 0.25 * tol1_score + 0.25 * tol2_score + 0.15 * len_score + 0.15 * mc_score + 0.20 * fdr_score
