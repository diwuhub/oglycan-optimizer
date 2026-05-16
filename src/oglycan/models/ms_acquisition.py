"""MS acquisition quality sub-model.

Scores resolution, AGC, injection time, and dynamic exclusion with
explicit duty-cycle trade-offs.

Key trade-off: higher resolution = fewer scans/sec. Longer injection
= better S/N but fewer spectra per cycle. These create a scan_rate
factor that penalizes greedy settings.

Literature basis:
  - Kelstrup et al., J Proteome Res 2012: Orbitrap resolution vs duty cycle
  - Thermo Fisher Orbitrap Exploris documentation: transient time scaling
"""

import math

from ._math import gauss as _gauss


def score(ms: dict) -> float:
    """Score MS acquisition quality.

    Parameters (from ms dict):
        resolution_ms1, resolution_ms2: Orbitrap resolution settings
        agc_target_ms1, agc_target_ms2: automatic gain control targets
        max_injection_ms1, max_injection_ms2: max injection times in ms
        dynamic_exclusion_s: dynamic exclusion window in seconds

    Returns:
        float in [0, 1]. Higher = better acquisition quality.
    """
    res1 = ms.get("resolution_ms1", 120000)
    res2 = ms.get("resolution_ms2", 30000)
    agc1 = ms.get("agc_target_ms1", 1e6)
    agc2 = ms.get("agc_target_ms2", 5e4)
    inj1 = ms.get("max_injection_ms1", 50)
    inj2 = ms.get("max_injection_ms2", 200)
    dyn_ex = ms.get("dynamic_exclusion_s", 30)

    # Resolution: diminishing returns past optimal
    res1_qual = _gauss(res1, 120000, 40000)
    res2_qual = _gauss(res2, 30000, 8000)

    # Scan rate penalty: Orbitrap transient time scales with resolution
    ms1_transient = res1 / 4000.0
    ms2_transient = res2 / 4000.0
    ms1_scan_ms = ms1_transient + inj1
    ms2_scan_ms = ms2_transient + inj2
    cycle_time_s = (ms1_scan_ms + 10 * ms2_scan_ms) / 1000.0
    scan_rate_score = _gauss(cycle_time_s, 2.5, 1.0)

    # AGC targets
    agc1_score = _gauss(math.log10(max(1, agc1)), 6.0, 0.5)
    agc2_score = _gauss(math.log10(max(1, agc2)), 4.7, 0.4)

    # Dynamic exclusion: 20-40s optimal
    dex_score = _gauss(dyn_ex, 30.0, 12.0)

    return (0.12 * res1_qual + 0.12 * res2_qual + 0.30 * scan_rate_score
            + 0.12 * agc1_score + 0.12 * agc2_score + 0.12 * dex_score
            + 0.10 * (res2_qual * scan_rate_score))
