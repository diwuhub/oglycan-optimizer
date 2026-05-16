"""Tests for canonical and MSFragger pilot ingestion."""

import json

import pytest

from oglycan.pilot import (
    load_pilot_byonic,
    load_pilot_canonical,
    load_pilot_msfragger,
    load_pilot_o_pair,
    validate_pilot,
)


def test_load_pilot_canonical_round_trip(tmp_path):
    payload = {
        "source": "canonical",
        "glycoprotein_name": "Etanercept",
        "sites": [
            {
                "pos": 266,
                "aa": "T",
                "observed_localization": 0.3,
                "n_spectra": 2,
                "observed_glycoforms": ["Core1_Sia", "Core2"],
            }
        ],
        "metadata": {"path": "pilot.json"},
    }
    pilot_path = tmp_path / "pilot.json"
    pilot_path.write_text(json.dumps(payload), encoding="utf-8")

    assert load_pilot_canonical(str(pilot_path)) == payload


def test_validate_pilot_accepts_valid_result_and_flags_common_errors():
    valid = {
        "source": "canonical",
        "glycoprotein_name": "Etanercept",
        "sites": [
            {
                "pos": 266,
                "aa": "T",
                "observed_localization": 0.95,
                "n_spectra": 3,
                "observed_glycoforms": ["Core1_Sia"],
            }
        ],
        "metadata": {"run_date": "2026-04-17"},
    }
    assert validate_pilot(valid) == []

    invalid = {
        "glycoprotein_name": "Etanercept",
        "sites": [
            {
                "aa": "T",
                "observed_localization": 0.95,
                "n_spectra": 3,
                "observed_glycoforms": ["Core1_Sia"],
            }
        ],
        "metadata": {},
    }
    errors = validate_pilot(invalid)
    assert errors
    assert any("missing source" in error for error in errors)
    assert any("missing 'pos'" in error for error in errors)


def test_load_pilot_msfragger_parses_minimal_tsv(tmp_path):
    tsv = "\n".join(
        [
            "\t".join(
                [
                    "Peptide",
                    "Observed Glycan Mass",
                    "Total Glycan Composition",
                    "Glycan Score",
                    "Best Localization",
                    "Protein Start",
                    "Protein End",
                    "Observed Modifications",
                ]
            ),
            "\t".join(
                [
                    "TDFT[+203.079]MPD",
                    "203.079",
                    "HexNAc(1)Hex(1)",
                    "15",
                    "0.8",
                    "100",
                    "106",
                    "4",
                ]
            ),
            "\t".join(
                [
                    "TDFT[+203.079]MPD",
                    "203.079",
                    "HexNAc(1)Hex(1)NeuAc(1)",
                    "25",
                    "0.9",
                    "100",
                    "106",
                    "4",
                ]
            ),
        ]
    )
    pilot_path = tmp_path / "pilot.tsv"
    pilot_path.write_text(tsv, encoding="utf-8")

    result = load_pilot_msfragger(str(pilot_path), glycoprotein_name="Etanercept")

    assert result["source"] == "msfragger_glyco"
    assert result["glycoprotein_name"] == "Etanercept"
    assert result["metadata"]["rows_seen"] == 2
    assert result["sites"] == [
        {
            "pos": 103,
            "aa": "T",
            "observed_localization": 0.9,
            "n_spectra": 2,
            "observed_glycoforms": ["Core1", "Core1_Sia"],
        }
    ]


def test_load_pilot_byonic_stub_raises():
    with pytest.raises(NotImplementedError, match="Byonic parser pending"):
        load_pilot_byonic("byonic.csv")


def test_load_pilot_o_pair_happy_path(tmp_path):
    tsv = "\n".join(
        [
            "\t".join(
                [
                    "Protein Accession",
                    "Base Sequence",
                    "Site Specific Localization Probability",
                    "Total Glycan Composition",
                    "Start and End Residues In Protein",
                    "Glycan Site",
                ]
            ),
            "\t".join(
                [
                    "sp|ETANERCEPT|Etanercept",
                    "AASTP",
                    "0.82",
                    "HexNAc(1)Hex(1)",
                    "[234-238]",
                    "4",
                ]
            ),
            "\t".join(
                [
                    "sp|ETANERCEPT|Etanercept",
                    "AASTP",
                    "0.91",
                    "HexNAc(1)Hex(1)NeuAc(1)",
                    "[234-238]",
                    "4",
                ]
            ),
            "\t".join(
                [
                    "sp|ETANERCEPT|Etanercept",
                    "TASGP",
                    "0.75",
                    "HexNAc(2)Hex(1)",
                    "[239-243]",
                    "3",
                ]
            ),
            "\t".join(
                [
                    "sp|ETANERCEPT|Etanercept",
                    "AATTK",
                    "0.88",
                    "HexNAc(1)",
                    "[263-267]",
                    "4",
                ]
            ),
        ]
    )
    pilot_path = tmp_path / "opair.tsv"
    pilot_path.write_text(tsv, encoding="utf-8")

    result = load_pilot_o_pair(str(pilot_path), glycoprotein_name="Etanercept")

    assert result["source"] == "o_pair"
    assert result["glycoprotein_name"] == "Etanercept"
    assert result["sites"] == [
        {
            "pos": 237,
            "aa": "T",
            "observed_localization": 0.91,
            "n_spectra": 2,
            "observed_glycoforms": ["Core1", "Core1_Sia"],
        },
        {
            "pos": 241,
            "aa": "S",
            "observed_localization": 0.75,
            "n_spectra": 1,
            "observed_glycoforms": ["Core2"],
        },
        {
            "pos": 266,
            "aa": "T",
            "observed_localization": 0.88,
            "n_spectra": 1,
            "observed_glycoforms": ["Core1_GalNAc"],
        },
    ]


def test_load_pilot_o_pair_missing_column_raises_value_error(tmp_path):
    tsv = "\n".join(
        [
            "\t".join(
                [
                    "Protein Accession",
                    "Base Sequence",
                    "Total Glycan Composition",
                    "Start and End Residues In Protein",
                    "Glycan Site",
                ]
            ),
            "\t".join(
                [
                    "Etanercept",
                    "AASTP",
                    "HexNAc(1)Hex(1)",
                    "[234-238]",
                    "4",
                ]
            ),
        ]
    )
    pilot_path = tmp_path / "opair_missing.tsv"
    pilot_path.write_text(tsv, encoding="utf-8")

    with pytest.raises(ValueError) as excinfo:
        load_pilot_o_pair(str(pilot_path), glycoprotein_name="Etanercept")

    message = str(excinfo.value)
    assert "Missing required O-Pair columns" in message
    assert "Site Specific Localization Probability" in message
    assert "Site Localization Probability" in message
    assert "Localization Probability" in message


def test_load_pilot_o_pair_aggregates_multi_row_sites(tmp_path):
    tsv = "\n".join(
        [
            "\t".join(
                [
                    "Organism",
                    "Base Sequence",
                    "Localization Probability",
                    "Glycan Composition",
                    "Start Residue In Protein",
                    "Localized Glycans",
                ]
            ),
            "\t".join(
                [
                    "Etanercept",
                    "TASGP",
                    "0.42",
                    "HexNAc(1)Hex(1)",
                    "239",
                    "S3:HexNAc(1)Hex(1)",
                ]
            ),
            "\t".join(
                [
                    "Etanercept",
                    "TASGP",
                    "0.87",
                    "HexNAc(1)Hex(1)NeuAc(1)",
                    "239",
                    "S3:HexNAc(1)Hex(1)NeuAc(1)",
                ]
            ),
        ]
    )
    pilot_path = tmp_path / "opair_aggregate.tsv"
    pilot_path.write_text(tsv, encoding="utf-8")

    result = load_pilot_o_pair(str(pilot_path), glycoprotein_name="Etanercept")

    assert result["sites"] == [
        {
            "pos": 241,
            "aa": "S",
            "observed_localization": 0.87,
            "n_spectra": 2,
            "observed_glycoforms": ["Core1", "Core1_Sia"],
        }
    ]
