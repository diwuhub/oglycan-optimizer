"""Tests for Phase 4 prediction scaffolding."""

import io
import urllib.error
import urllib.parse
from unittest.mock import patch

import pytest

from oglycan.predict import (
    _parse_netoglyc_response,
    parse_fasta,
    predict_netoglyc,
    scan_st_positions,
    suggest_catalog,
    validate_predicted_sites,
)

VALID_NETOGLYC_RESPONSE = "\n".join(
    [
        "# NetOGlyc 4.0",
        "# columns: position aa score binary_prediction",
        "# query: toy_protein",
        "1\tM\t0.01\t-",
        "2\tS\t0.91\t+",
        "3\tT\t0.62\t+",
        "4\tP\t0.02\t-",
        "5\tM\t0.03\t-",
        "6\tS\t0.40\t-",
        "7\tA\t0.01\t-",
        "8\tG\t0.00\t-",
    ]
)

THRESHOLD_NETOGLYC_RESPONSE = "\n".join(
    [
        "# NetOGlyc 4.0",
        "# columns: position aa score binary_prediction",
        "# query: toy_protein",
        "1\tM\t0.01\t-",
        "2\tS\t0.90\t+",
        "3\tT\t0.30\t-",
        "4\tP\t0.02\t-",
        "5\tM\t0.03\t-",
        "6\tS\t0.60\t+",
    ]
)


class _MockHTTPResponse:
    def __init__(self, body: str, *, status: int = 200):
        self._body = io.BytesIO(body.encode("utf-8"))
        self.status = status

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        return False

    def getcode(self):
        return self.status

    def read(self):
        return self._body.read()


class _AssertingNetOGlycOpener:
    def __init__(self, response_text: str, *, expected_sequence: str):
        self.expected_sequence = expected_sequence
        self.response_text = response_text
        self.calls: list[tuple[object, float]] = []

    def open(self, request, timeout=120.0):
        self.calls.append((request, timeout))
        assert request.full_url == "https://services.healthtech.dtu.dk/services/NetOGlyc-4.0/"
        assert request.get_method() == "POST"
        form = urllib.parse.parse_qs(request.data.decode("utf-8"), strict_parsing=True)
        assert form == {
            "SEQPASTE": [self.expected_sequence],
            "outputformat": ["short"],
        }
        return _MockHTTPResponse(self.response_text)


def test_parse_fasta_single_record_normalizes_uppercase():
    fasta_text = ">toy\nmstp\n"
    assert parse_fasta(fasta_text) == [("toy", "MSTP")]


def test_parse_fasta_multi_record_handles_multiline_sequences():
    fasta_text = ">first\nms\n tp\n>second\naa\nst\n"
    assert parse_fasta(fasta_text) == [("first", "MSTP"), ("second", "AAST")]


def test_parse_fasta_empty_input_returns_no_records():
    assert parse_fasta("") == []


def test_scan_st_positions_returns_expected_candidates():
    assert scan_st_positions("MSTP") == [
        {"pos": 2, "aa": "S", "p_glycosite": None, "source": "scan_st"},
        {"pos": 3, "aa": "T", "p_glycosite": None, "source": "scan_st"},
    ]


def test_scan_st_positions_covers_etanercept_stalk_sites():
    stalk_subsequence = "TASATATASAATASAATATAASATAASAAT"
    predicted_sites = scan_st_positions(stalk_subsequence)
    assert [site["pos"] + 236 for site in predicted_sites] == [
        237,
        239,
        241,
        243,
        245,
        248,
        250,
        253,
        255,
        258,
        260,
        263,
        266,
    ]
    assert [site["aa"] for site in predicted_sites] == [
        "T",
        "S",
        "T",
        "T",
        "S",
        "T",
        "S",
        "T",
        "T",
        "S",
        "T",
        "S",
        "T",
    ]


def test_suggest_catalog_scan_st_returns_empty_sites_and_candidates():
    catalog = suggest_catalog("MSTPAAAPSTPAAA", "toy_protein")
    assert catalog["glycoprotein"]["name"] == "toy_protein"
    assert catalog["sites"] == []
    assert len(catalog["predicted_sites"]) == 4


def test_suggest_catalog_netoglyc_passes_threshold_to_predictor():
    predicted_sites = [{"pos": 2, "aa": "S", "p_glycosite": 0.7, "source": "netoglyc_4.0"}]

    with patch("oglycan.predict.predict_netoglyc", return_value=predicted_sites) as mock_predict:
        catalog = suggest_catalog(
            "MSTP",
            "toy_protein",
            predictor="netoglyc",
            threshold=0.7,
        )

    mock_predict.assert_called_once_with("MSTP", threshold=0.7)
    assert catalog["predicted_sites"] == predicted_sites


def test_validate_predicted_sites_accepts_valid_entries():
    assert validate_predicted_sites(
        [
            {"pos": 2, "aa": "S", "p_glycosite": None, "source": "scan_st"},
            {"pos": 3, "aa": "T", "p_glycosite": 0.42, "source": "netoglyc_4.0"},
        ]
    ) == []


def test_validate_predicted_sites_catches_missing_field():
    errors = validate_predicted_sites(
        [{"pos": 2, "aa": "S", "p_glycosite": None}]
    )
    assert any("missing field" in error and "'source'" in error for error in errors)


def test_validate_predicted_sites_catches_bad_aa_and_probability():
    errors = validate_predicted_sites(
        [
            {"pos": 2, "aa": "A", "p_glycosite": None, "source": "scan_st"},
            {"pos": 3, "aa": "T", "p_glycosite": 1.2, "source": "scan_st"},
        ]
    )
    assert any("aa must be one of" in error for error in errors)
    assert any("out of [0, 1]" in error for error in errors)


def test_netoglyc_parses_valid_tabular_response():
    assert _parse_netoglyc_response(VALID_NETOGLYC_RESPONSE) == [
        {"pos": 2, "aa": "S", "p_glycosite": 0.91, "source": "netoglyc_4.0"},
        {"pos": 3, "aa": "T", "p_glycosite": 0.62, "source": "netoglyc_4.0"},
        {"pos": 6, "aa": "S", "p_glycosite": 0.4, "source": "netoglyc_4.0"},
    ]


def test_netoglyc_raises_value_error_on_malformed_line():
    malformed = "\n".join(
        [
            "# NetOGlyc 4.0",
            "# columns: position aa score binary_prediction",
            "# query: toy_protein",
            "4\tS",
        ]
    )

    with pytest.raises(ValueError, match="tab-separated"):
        _parse_netoglyc_response(malformed)

    with pytest.raises(ValueError) as exc_info:
        _parse_netoglyc_response(malformed)

    assert "4\tS" in str(exc_info.value)


def test_netoglyc_raises_when_no_data_rows():
    opener = _AssertingNetOGlycOpener(
        "\n".join(
            [
                "# NetOGlyc 4.0",
                "# columns: position aa score binary_prediction",
                "# query: toy_protein",
            ]
        ),
        expected_sequence="MSTP",
    )

    with pytest.raises(RuntimeError, match="no predictions"):
        predict_netoglyc("MSTP", opener=opener)


def test_netoglyc_http_404_raises_runtime_error():
    http_error = urllib.error.HTTPError(
        "https://services.healthtech.dtu.dk/services/NetOGlyc-4.0/",
        404,
        "Not Found",
        hdrs=None,
        fp=io.BytesIO(b"missing"),
    )

    with patch("oglycan.predict.urllib.request.urlopen", side_effect=http_error):
        with pytest.raises(RuntimeError, match="HTTP 404"):
            predict_netoglyc("MSTP")


def test_netoglyc_happy_path_end_to_end_with_mock_opener():
    response_text = (
        VALID_NETOGLYC_RESPONSE
        .replace("3\tT\t0.62\t+", "3\tT\t0.72\t+")
        .replace("6\tS\t0.40\t-", "6\tS\t0.60\t+")
    )
    opener = _AssertingNetOGlycOpener(response_text, expected_sequence="MSTPMS")

    predicted_sites = predict_netoglyc("MSTPMS", opener=opener)

    assert opener.calls
    assert predicted_sites == [
        {"pos": 2, "aa": "S", "p_glycosite": 0.91, "source": "netoglyc_4.0"},
        {"pos": 3, "aa": "T", "p_glycosite": 0.72, "source": "netoglyc_4.0"},
        {"pos": 6, "aa": "S", "p_glycosite": 0.6, "source": "netoglyc_4.0"},
    ]


def test_netoglyc_respects_threshold():
    opener = _AssertingNetOGlycOpener(THRESHOLD_NETOGLYC_RESPONSE, expected_sequence="MSTPMS")

    predicted_sites = predict_netoglyc("MSTPMS", threshold=0.5, opener=opener)

    assert predicted_sites == [
        {"pos": 2, "aa": "S", "p_glycosite": 0.9, "source": "netoglyc_4.0"},
        {"pos": 6, "aa": "S", "p_glycosite": 0.6, "source": "netoglyc_4.0"},
    ]
