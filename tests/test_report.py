"""Tests for JSON report audit metadata."""

import hashlib
import json

from oglycan.method_card import _read_package_version
from oglycan.report import json_report


def test_json_report_includes_audit_fields(tmp_path):
    catalog_path = tmp_path / "catalog.json"
    params_path = tmp_path / "params.json"
    catalog_path.write_text('{"catalog":"known"}', encoding="utf-8")
    params_path.write_text('{"params":"known"}', encoding="utf-8")

    result = {
        "molecule": "Example",
        "total_sites": 0,
        "sites_localized": 0,
        "localization_threshold": 0.75,
        "composite_score": 0.5,
        "sub_model_scores": {},
        "composite_breakdown": {},
        "site_results": [],
    }
    full_params = {"fragmentation": {"mode": "EThcD"}}

    payload = json.loads(
        json_report(
            result,
            source_catalog=str(catalog_path),
            source_params=str(params_path),
            full_params=full_params,
        )
    )

    assert payload["package_version"] == _read_package_version()
    assert payload["source_catalog_path"] == str(catalog_path)
    assert payload["source_catalog_sha256"] == hashlib.sha256(catalog_path.read_bytes()).hexdigest()
    assert payload["source_params_path"] == str(params_path)
    assert payload["source_params_sha256"] == hashlib.sha256(params_path.read_bytes()).hexdigest()
    assert payload["acquisition_params"] == full_params
    assert payload["generated_at"].endswith("Z")
