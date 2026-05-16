"""Pilot-data ingestion helpers.

Canonical JSON is the universal ingestion path for pilot evidence consumed by
the optimizer. This module also supports narrow MSFragger-Glyco / FragPipe and
O-Pair / MetaMorpheus TSV parsers. Byonic remains intentionally unimplemented;
use canonical JSON input when its output needs to be ingested.
"""

from __future__ import annotations

import csv
import json
import math
import os
import re

# Canonical per-site pilot evidence.
# - pos: 1-indexed protein position
# - aa:  "S" or "T"
# - observed_localization: float in [0,1]; best/max site-localization
#        probability from the pilot run (EThcD -> c/z-ion assignment
#        confidence). NaN-free; use None if site was not sampled.
# - n_spectra: int; number of PSMs that localized to this site (0 if not seen)
# - observed_glycoforms: list[str]; canonical core_types actually observed
PilotSite = dict

GLYCAN_COMPOSITION_TO_CORE = {
    "HexNAc(1)": "Core1_GalNAc",
    "HexNAc(1)Hex(1)": "Core1",
    "HexNAc(1)Hex(1)NeuAc(1)": "Core1_Sia",
    "HexNAc(1)Hex(1)NeuAc(2)": "Core1_Sia",
    "HexNAc(2)Hex(1)": "Core2",
    "HexNAc(2)Hex(2)": "Core2",
    "HexNAc(2)Hex(2)NeuAc(1)": "Core2",
    "HexNAc(2)Hex(2)Fuc(1)": "Core2_Fuc",
    "HexNAc(2)Hex(2)NeuAc(1)Fuc(1)": "Core2_Fuc",
}

_VALID_SOURCES = {"msfragger_glyco", "byonic", "o_pair", "canonical"}
_MSFRAGGER_COLUMN_GROUPS = {
    "peptide": ["Peptide"],
    "glycan_mass": ["Observed Glycan Mass", "Glycan Mass"],
    "glycan_composition": ["Total Glycan Composition"],
    "glycan_score": ["Glycan Score", "Protein Glycan Score"],
    "localization": ["Best Localization", "Localization Probability"],
    "protein_start": ["Protein Start"],
    "protein_end": ["Protein End"],
    "site_offset": [
        "Observed Modifications",
        "Assigned Modifications",
        "Site Position in Peptide",
        "Glycan Site",
    ],
}
_OPAIR_COLUMN_GROUPS = {
    "protein_name": ["Protein Accession", "Organism"],
    "peptide": ["Base Sequence", "Full Sequence", "Peptide"],
    "localization": [
        "Site Specific Localization Probability",
        "Site Localization Probability",
        "Localization Probability",
    ],
    "glycan_composition": [
        "Total Glycan Composition",
        "Glycan Composition",
        "Oxonium Ion Composition",
    ],
    "protein_start": [
        "Start and End Residues In Protein",
        "Start Residue In Protein",
    ],
    "site_offset": ["Glycan Site", "Localized Glycans"],
}


def load_pilot_canonical(path: str) -> dict:
    """Load canonical pilot JSON without mutating its contents."""
    with open(path) as f:
        return json.load(f)


def validate_pilot(result: dict) -> list[str]:
    """Validate canonical pilot-result structure."""
    errors = []
    if not isinstance(result, dict):
        return ["pilot result must be a dict"]

    source = result.get("source")
    if not source:
        errors.append("missing source")
    elif source not in _VALID_SOURCES:
        errors.append(f"invalid source: {source!r}")

    glycoprotein_name = result.get("glycoprotein_name")
    if not isinstance(glycoprotein_name, str) or not glycoprotein_name.strip():
        errors.append("missing glycoprotein_name")

    sites = result.get("sites")
    if not isinstance(sites, list):
        errors.append("missing sites")
        sites = []

    metadata = result.get("metadata")
    if metadata is None:
        errors.append("missing metadata")
    elif not isinstance(metadata, dict):
        errors.append("metadata must be a dict")

    for idx, site in enumerate(sites):
        prefix = f"site[{idx}]"
        if not isinstance(site, dict):
            errors.append(f"{prefix}: must be a dict")
            continue

        pos = site.get("pos")
        if pos is None:
            errors.append(f"{prefix}: missing 'pos'")
        elif not isinstance(pos, int) or pos < 1:
            errors.append(f"{prefix}: invalid 'pos' {pos!r}")

        aa = site.get("aa")
        if aa is None:
            errors.append(f"{prefix}: missing 'aa'")
        elif aa not in {"S", "T"}:
            errors.append(f"{prefix}: invalid 'aa' {aa!r}")

        if "observed_localization" not in site:
            errors.append(f"{prefix}: missing 'observed_localization'")
        else:
            loc = site["observed_localization"]
            if loc is not None:
                if not isinstance(loc, (int, float)):
                    errors.append(f"{prefix}: observed_localization must be numeric or None")
                elif math.isnan(float(loc)):
                    errors.append(f"{prefix}: observed_localization must not be NaN")
                elif not (0.0 <= float(loc) <= 1.0):
                    errors.append(f"{prefix}: observed_localization {loc!r} out of [0, 1]")

        if "n_spectra" not in site:
            errors.append(f"{prefix}: missing 'n_spectra'")
        else:
            n_spectra = site["n_spectra"]
            if not isinstance(n_spectra, int) or n_spectra < 0:
                errors.append(f"{prefix}: invalid 'n_spectra' {n_spectra!r}")

        glycoforms = site.get("observed_glycoforms")
        if glycoforms is None:
            errors.append(f"{prefix}: missing 'observed_glycoforms'")
        elif not isinstance(glycoforms, list):
            errors.append(f"{prefix}: observed_glycoforms must be a list")
        elif any(not isinstance(item, str) for item in glycoforms):
            errors.append(f"{prefix}: observed_glycoforms entries must be strings")

    return errors


def load_pilot_msfragger(
    path: str,
    glycoprotein_name: str,
    sequence: str | None = None,
) -> dict:
    """Parse a narrow MSFragger-Glyco / FragPipe psm.tsv into canonical form."""
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError("MSFragger TSV has no header row")

        columns, missing = _resolve_msfragger_columns(reader.fieldnames)
        if missing:
            raise ValueError(
                "Missing required MSFragger columns: " + ", ".join(missing)
            )

        grouped = {}
        rows_seen = 0

        for row_index, row in enumerate(reader, start=2):
            if not any((value or "").strip() for value in row.values()):
                continue

            rows_seen += 1
            peptide = (row.get(columns["peptide"]) or "").strip()
            if not peptide:
                raise ValueError(f"row {row_index}: missing peptide sequence")

            residues = _peptide_residues(peptide)
            if not residues:
                raise ValueError(f"row {row_index}: could not parse peptide residues")

            protein_start = _parse_int(
                row.get(columns["protein_start"]), "Protein Start", row_index
            )
            protein_end = _parse_int(
                row.get(columns["protein_end"]), "Protein End", row_index
            )
            offset = _extract_site_offset(
                row.get(columns["site_offset"]), peptide, row_index
            )
            if not (1 <= offset <= len(residues)):
                raise ValueError(
                    f"row {row_index}: site offset {offset} outside peptide length {len(residues)}"
                )

            pos = protein_start + offset - 1
            if pos > protein_end:
                raise ValueError(
                    f"row {row_index}: site position {pos} exceeds Protein End {protein_end}"
                )

            if sequence is not None:
                if pos > len(sequence):
                    raise ValueError(
                        f"row {row_index}: site position {pos} outside provided sequence length"
                    )
                aa = sequence[pos - 1].upper()
            else:
                aa = residues[offset - 1]

            if aa not in {"S", "T"}:
                raise ValueError(
                    f"row {row_index}: localized amino acid must be S/T, got {aa!r}"
                )

            observed_localization = _parse_optional_float(
                row.get(columns["localization"]), row_index, columns["localization"]
            )
            composition = _normalize_composition(row.get(columns["glycan_composition"]))
            core_type = GLYCAN_COMPOSITION_TO_CORE.get(composition, "Unknown")

            key = (pos, aa)
            site = grouped.setdefault(
                key,
                {
                    "pos": pos,
                    "aa": aa,
                    "observed_localization": None,
                    "n_spectra": 0,
                    "observed_glycoforms": set(),
                },
            )
            site["n_spectra"] += 1
            if observed_localization is not None:
                current = site["observed_localization"]
                site["observed_localization"] = (
                    observed_localization
                    if current is None
                    else max(current, observed_localization)
                )
            site["observed_glycoforms"].add(core_type)

    result = {
        "source": "msfragger_glyco",
        "glycoprotein_name": glycoprotein_name,
        "sites": [
            {
                **site,
                "observed_glycoforms": sorted(site["observed_glycoforms"]),
            }
            for _, site in sorted(grouped.items())
        ],
        "metadata": {
            "path": os.path.abspath(path),
            "parser": "msfragger_psm_tsv",
            "rows_seen": rows_seen,
        },
    }
    errors = validate_pilot(result)
    if errors:
        raise ValueError("Invalid canonical pilot result: " + "; ".join(errors))
    return result


def load_pilot_byonic(*_args, **_kwargs):
    raise NotImplementedError(
        "Byonic parser pending — use canonical JSON input or MSFragger TSV"
    )
def load_pilot_o_pair(
    path: str,
    glycoprotein_name: str,
    sequence: str | None = None,
) -> dict:
    """
    Parse O-Pair (MetaMorpheus) per-PSM TSV output into canonical PilotResult.
    """
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError("O-Pair TSV has no header row")

        columns, missing = _resolve_o_pair_columns(reader.fieldnames)
        if missing:
            raise ValueError("Missing required O-Pair columns: " + ", ".join(missing))

        grouped = {}
        rows_seen = 0

        for row_index, row in enumerate(reader, start=2):
            if not any((value or "").strip() for value in row.values()):
                continue

            rows_seen += 1
            protein_name = (row.get(columns["protein_name"]) or "").strip()
            if not _matches_glycoprotein_name(protein_name, glycoprotein_name):
                continue

            peptide = (row.get(columns["peptide"]) or "").strip()
            if not peptide:
                raise ValueError(f"row {row_index}: missing peptide sequence")

            residues = _peptide_residues(peptide)
            if not residues:
                raise ValueError(f"row {row_index}: could not parse peptide residues")

            protein_start = _parse_o_pair_protein_start(
                row.get(columns["protein_start"]),
                row_index,
                columns["protein_start"],
            )
            offset = _extract_o_pair_site_offset(
                row.get(columns["site_offset"]),
                peptide,
                row_index,
                columns["site_offset"],
            )
            if not (1 <= offset <= len(residues)):
                raise ValueError(
                    f"row {row_index}: site offset {offset} outside peptide length {len(residues)}"
                )

            pos = protein_start + offset - 1
            if sequence is not None:
                if pos > len(sequence):
                    raise ValueError(
                        f"row {row_index}: site position {pos} outside provided sequence length"
                    )
                aa = sequence[pos - 1].upper()
            else:
                aa = residues[offset - 1]

            if aa not in {"S", "T"}:
                raise ValueError(
                    f"row {row_index}: localized amino acid must be S/T, got {aa!r}"
                )

            observed_localization = _parse_optional_float(
                row.get(columns["localization"]), row_index, columns["localization"]
            )
            composition = _normalize_composition(row.get(columns["glycan_composition"]))
            core_type = GLYCAN_COMPOSITION_TO_CORE.get(composition, "Unknown")

            key = (pos, aa)
            site = grouped.setdefault(
                key,
                {
                    "pos": pos,
                    "aa": aa,
                    "observed_localization": None,
                    "n_spectra": 0,
                    "observed_glycoforms": set(),
                },
            )
            site["n_spectra"] += 1
            if observed_localization is not None:
                current = site["observed_localization"]
                site["observed_localization"] = (
                    observed_localization
                    if current is None
                    else max(current, observed_localization)
                )
            site["observed_glycoforms"].add(core_type)

    result = {
        "source": "o_pair",
        "glycoprotein_name": glycoprotein_name,
        "sites": [
            {
                **site,
                "observed_glycoforms": sorted(site["observed_glycoforms"]),
            }
            for _, site in sorted(grouped.items())
        ],
        "metadata": {
            "path": os.path.abspath(path),
            "parser": "o_pair_psm_tsv",
            "rows_seen": rows_seen,
        },
    }
    errors = validate_pilot(result)
    if errors:
        raise ValueError("Invalid canonical pilot result: " + "; ".join(errors))
    return result


def _resolve_msfragger_columns(fieldnames: list[str]) -> tuple[dict, list[str]]:
    return _resolve_columns(fieldnames, _MSFRAGGER_COLUMN_GROUPS)


def _resolve_o_pair_columns(fieldnames: list[str]) -> tuple[dict, list[str]]:
    return _resolve_columns(fieldnames, _OPAIR_COLUMN_GROUPS)


def _resolve_columns(
    fieldnames: list[str],
    column_groups: dict[str, list[str]],
) -> tuple[dict, list[str]]:
    columns = {}
    missing = []
    available = set(fieldnames)
    for label, candidates in column_groups.items():
        column = next((candidate for candidate in candidates if candidate in available), None)
        columns[label] = column
        if column is None:
            missing.append("/".join(candidates))
    return columns, missing


def _parse_int(value: str | None, field_name: str, row_index: int) -> int:
    text = (value or "").strip()
    if not text:
        raise ValueError(f"row {row_index}: missing {field_name}")
    try:
        return int(text)
    except ValueError as exc:
        raise ValueError(f"row {row_index}: invalid {field_name} {text!r}") from exc


def _parse_optional_float(value: str | None, row_index: int, field_name: str) -> float | None:
    text = (value or "").strip()
    if not text:
        return None
    try:
        parsed = float(text)
    except ValueError as exc:
        raise ValueError(f"row {row_index}: invalid {field_name} {text!r}") from exc
    if math.isnan(parsed):
        return None
    if not (0.0 <= parsed <= 1.0):
        raise ValueError(f"row {row_index}: {field_name} {parsed!r} out of [0, 1]")
    return parsed


def _normalize_composition(value: str | None) -> str:
    return "".join((value or "").split())


def _matches_glycoprotein_name(observed: str | None, expected: str) -> bool:
    observed_text = (observed or "").strip().lower()
    expected_text = expected.strip().lower()
    return bool(observed_text and expected_text and expected_text in observed_text)


def _parse_o_pair_protein_start(
    value: str | None, row_index: int, field_name: str
) -> int:
    text = (value or "").strip()
    if not text:
        raise ValueError(f"row {row_index}: missing {field_name}")
    if field_name == "Start and End Residues In Protein":
        match = re.search(r"\[\s*(\d+)", text)
        if not match:
            raise ValueError(f"row {row_index}: invalid {field_name} {text!r}")
        return int(match.group(1))
    return _parse_int(value, field_name, row_index)


def _extract_o_pair_site_offset(
    value: str | None,
    peptide: str,
    row_index: int,
    field_name: str,
) -> int:
    text = (value or "").strip()
    if not text:
        raise ValueError(f"row {row_index}: missing {field_name}")
    if field_name == "Glycan Site":
        return _parse_int(value, field_name, row_index)

    candidates = []
    for pattern in (
        r"\b[A-Z](\d+)\b",
        r"\b[A-Z]\((\d+)\)",
        r"\b(\d+)(?=:[A-Z]\b)",
        r"\bsite\D*(\d+)\b",
        r"\bposition\D*(\d+)\b",
    ):
        candidates.extend(int(match) for match in re.findall(pattern, text, flags=re.IGNORECASE))

    if candidates:
        residues = _peptide_residues(peptide)
        for candidate in candidates:
            if 1 <= candidate <= len(residues):
                return candidate
        return candidates[0]

    return _extract_site_offset(value, peptide, row_index)


def _extract_site_offset(value: str | None, peptide: str, row_index: int) -> int:
    text = (value or "").strip()
    residues = _peptide_residues(peptide)

    if text:
        candidates = [int(match) for match in re.findall(r"\d+", text)]
        for candidate in candidates:
            if 1 <= candidate <= len(residues):
                return candidate
        if candidates:
            return candidates[0]

    annotated = _annotated_site_offsets(peptide)
    if len(annotated) == 1:
        return annotated[0]
    if not annotated:
        raise ValueError(f"row {row_index}: could not determine localized site offset")
    raise ValueError(
        f"row {row_index}: multiple annotated sites found without an explicit offset"
    )


def _annotated_site_offsets(peptide: str) -> list[int]:
    offsets = []
    residue_index = 0
    in_mod = False
    for idx, char in enumerate(peptide):
        if char == "[":
            in_mod = True
            continue
        if char == "]":
            in_mod = False
            continue
        if in_mod:
            continue
        if char.isalpha() and char == char.upper():
            residue_index += 1
            if idx + 1 < len(peptide) and peptide[idx + 1] == "[":
                offsets.append(residue_index)
    return offsets


def _peptide_residues(peptide: str) -> list[str]:
    residues = []
    in_mod = False
    for char in peptide:
        if char == "[":
            in_mod = True
            continue
        if char == "]":
            in_mod = False
            continue
        if in_mod:
            continue
        if char.isalpha() and char == char.upper():
            residues.append(char)
    return residues
