"""OpeRATOR enzyme preprocessing sub-model.

Scores enzyme preprocessing strategy for O-glycopeptide analysis.

Key insight: OpeRATOR + SialEXO is the gold-standard combination.
OpeRATOR cleaves at O-glycosites, providing direct site localization.
Without SialEXO pretreatment, OpeRATOR has greatly reduced activity
on sialylated glycoforms (+17% yield with SialEXO).

Literature basis:
  - Trastoy et al., PNAS 2020: OpeRATOR mechanism
  - Genovis application note AN-0042: SialEXO prerequisite for OpeRATOR
"""

from ._math import gauss as _gauss


def score(enzyme: dict) -> float:
    """Score enzyme preprocessing quality.

    Parameters (from enzyme dict):
        primary: enzyme name string (e.g. "OpeRATOR", "OglyZOR")
        use_sialexo: whether SialEXO pretreatment is used
        use_pngasef: whether PNGaseF is used to remove N-glycan interference
        digestion_time_hours: digestion duration (default 4)
        temperature_C: digestion temperature (default 37)

    Returns:
        float in [0, 0.35]. This is a bonus score, not a full 0-1 range.
    """
    primary = enzyme.get("primary", "").lower()
    use_sialexo = enzyme.get("use_sialexo", False)
    use_pngasef = enzyme.get("use_pngasef", False)
    dig_time = enzyme.get("digestion_time_hours", 4)
    dig_temp = enzyme.get("temperature_C", 37)

    bonus = 0.0
    if "operator" in primary:
        bonus = 0.25 if use_sialexo else 0.08
        time_factor = _gauss(dig_time, 4.0, 2.0)
        bonus *= (0.4 + 0.6 * time_factor)
        temp_factor = _gauss(dig_temp, 37.0, 4.0)
        bonus *= (0.3 + 0.7 * temp_factor)

    elif "oglyzor" in primary:
        bonus = 0.15 if use_sialexo else 0.05
        time_factor = _gauss(dig_time, 3.0, 1.5)
        bonus *= (0.4 + 0.6 * time_factor)

    if use_pngasef:
        bonus += 0.06

    return min(0.35, bonus)
