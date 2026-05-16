"""Composition layer: sub-models -> site scores -> composite.

Loads a site catalog (JSON), runs all six sub-models against acquisition
parameters, scores each site, and produces aggregate metrics.

This is the single entry point for evaluation. CLI and report modules
call core.evaluate().
"""

import json
import os

from .models import (
    etd_efficiency,
    hcd_fragments,
    lc_separation,
    ms_acquisition,
    operator_bonus,
    search_quality,
)
from .models._math import gauss as _gauss
from .site_scoring import score_site


def load_site_catalog(path: str) -> dict:
    """Load a glycoprotein site catalog from JSON."""
    with open(path) as f:
        return json.load(f)


def load_params(path: str) -> dict:
    """Load acquisition parameters from JSON."""
    with open(path) as f:
        return json.load(f)


def evaluate(site_catalog: dict, params: dict) -> dict:
    """Run full evaluation.

    Args:
        site_catalog: parsed site catalog (from sites/*.json).
        params: parsed acquisition parameters.

    Returns:
        dict with all scores, per-site results, and aggregate metrics.
    """
    enzyme = params.get("enzyme_preprocessing", {})
    frag = params.get("fragmentation", {})
    ms_acq = params.get("ms_acquisition", {})
    glycan_db = params.get("glycan_database", {})
    search = params.get("search_params", {})
    lc = params.get("lc_gradient", {})

    etd_eff = etd_efficiency.score(frag)
    hcd_frag = hcd_fragments.score(frag, glycan_db)
    op_bonus = operator_bonus.score(enzyme)
    ms_qual = ms_acquisition.score(ms_acq)
    lc_sep = lc_separation.score(lc)
    search_qual = search_quality.score(search)

    frag_mode = frag.get("mode", "EThcD").lower()
    if frag_mode != "ethcd":
        if frag_mode == "hcd":
            etd_eff *= 0.2
        elif frag_mode == "etd":
            hcd_frag *= 0.3
        else:
            etd_eff *= 0.5
            hcd_frag *= 0.5

    sites = site_catalog.get("sites", [])
    threshold = site_catalog.get("localization_threshold", 0.75)
    num_sites = len(sites)

    site_results = []
    for site in sites:
        sc = score_site(site, etd_eff, hcd_frag, op_bonus,
                        search_qual, lc_sep, ms_qual)
        site_results.append({
            "position": site["pos"],
            "amino_acid": site["aa"],
            "core_types": site["core_types"],
            "difficulty": site["difficulty"],
            "confidence": round(sc, 4),
            "pass": sc >= threshold,
        })

    site_scores = [s["confidence"] for s in site_results]
    sites_localized = sum(1 for s in site_results if s["pass"])
    localization_confidence = sum(site_scores) / max(1, len(site_scores))

    spectral_quality = (0.30 * etd_eff + 0.25 * hcd_frag
                        + 0.25 * ms_qual + 0.20 * lc_sep)

    db_cov = 0.0
    if glycan_db.get("include_core1"): db_cov += 0.30
    if glycan_db.get("include_core2"): db_cov += 0.25
    if glycan_db.get("include_sialylated"): db_cov += 0.20
    if glycan_db.get("include_fucosylated"): db_cov += 0.15
    if glycan_db.get("include_core3"): db_cov += 0.05
    max_sz = glycan_db.get("max_glycan_size", 4)
    sz_factor = _gauss(max_sz, 5, 1.5)
    glycan_diversity = min(1.0, hcd_frag * search_qual * db_cov
                           * (0.5 + 0.5 * sz_factor) * 2.5)

    bio_plaus = 0.40
    primary = enzyme.get("primary", "").lower()
    if "operator" in primary:
        bio_plaus += 0.20 if enzyme.get("use_sialexo") else 0.05
    elif "oglyzor" in primary:
        bio_plaus += 0.12 if enzyme.get("use_sialexo") else 0.03
    if enzyme.get("use_pngasef"):
        bio_plaus += 0.10
    if glycan_db.get("include_core1") and glycan_db.get("include_core2"):
        bio_plaus += 0.10
    if glycan_db.get("include_sialylated"):
        bio_plaus += 0.08
    if glycan_db.get("include_fucosylated"):
        bio_plaus += 0.07
    if glycan_db.get("include_fucosylated") and glycan_db.get("include_core2"):
        bio_plaus += 0.05
    bio_plaus = min(1.0, bio_plaus)

    sequence_coverage = sites_localized / max(1, num_sites)

    composite = (0.25 * spectral_quality + 0.25 * localization_confidence
                 + 0.20 * sequence_coverage + 0.15 * glycan_diversity
                 + 0.15 * bio_plaus)

    return {
        "molecule": site_catalog.get("glycoprotein", {}).get("name", "unknown"),
        "total_sites": num_sites,
        "sites_localized": sites_localized,
        "localization_threshold": threshold,
        "composite_score": round(composite, 4),
        "sub_model_scores": {
            "etd_efficiency": round(etd_eff, 4),
            "hcd_fragments": round(hcd_frag, 4),
            "operator_bonus": round(op_bonus, 4),
            "ms_acquisition": round(ms_qual, 4),
            "lc_separation": round(lc_sep, 4),
            "search_quality": round(search_qual, 4),
        },
        "composite_breakdown": {
            "spectral_quality": round(spectral_quality, 4),
            "localization_confidence": round(localization_confidence, 4),
            "sequence_coverage": round(sequence_coverage, 4),
            "glycan_diversity": round(glycan_diversity, 4),
            "bio_plausibility": round(bio_plaus, 4),
        },
        "site_results": site_results,
    }
