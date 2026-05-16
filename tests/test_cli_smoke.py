"""CLI smoke tests for top-level subcommands."""

import json
import os
import re
import subprocess
import sys


def repo_root():
    return os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def test_cli_tune_runs():
    result = subprocess.run(
        [sys.executable, "-m", "oglycan", "tune", "sites/etanercept.json"],
        capture_output=True,
        text=True,
        timeout=30,
        cwd=repo_root(),
    )
    assert result.returncode == 0
    assert "COMPOSITE SCORE" in result.stdout


def test_cli_tune_with_out(tmp_path):
    out_path = tmp_path / "t.json"
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "oglycan",
            "tune",
            "sites/etanercept.json",
            "--out",
            str(out_path),
        ],
        capture_output=True,
        text=True,
        timeout=30,
        cwd=repo_root(),
    )
    assert result.returncode == 0
    assert "COMPOSITE SCORE" in result.stdout
    payload = json.loads(out_path.read_text())
    assert "package_version" in payload


def test_cli_tune_json_mode():
    result = subprocess.run(
        [sys.executable, "-m", "oglycan", "tune", "sites/etanercept.json", "--json"],
        capture_output=True,
        text=True,
        timeout=30,
        cwd=repo_root(),
    )
    assert result.returncode == 0
    payload = json.loads(result.stdout)
    assert "composite_score" in payload


def test_cli_protein_shorthand():
    shorthand = subprocess.run(
        [sys.executable, "-m", "oglycan", "tune", "etanercept"],
        capture_output=True,
        text=True,
        timeout=30,
        cwd=repo_root(),
    )
    full_path = subprocess.run(
        [sys.executable, "-m", "oglycan", "tune", "sites/etanercept.json"],
        capture_output=True,
        text=True,
        timeout=30,
        cwd=repo_root(),
    )

    assert shorthand.returncode == 0
    assert full_path.returncode == 0
    assert shorthand.stdout == full_path.stdout
    assert shorthand.stderr == full_path.stderr


def test_cli_method_card_runs(tmp_path):
    out_path = tmp_path / "method_card.md"
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "oglycan",
            "method-card",
            "sites/etanercept.json",
            "--out",
            str(out_path),
        ],
        capture_output=True,
        text=True,
        timeout=30,
        cwd=repo_root(),
    )
    assert result.returncode == 0
    assert out_path.exists()
    assert out_path.read_text()


def test_cli_uncertainty_runs():
    result = subprocess.run(
        [sys.executable, "-m", "oglycan", "uncertainty", "sites/etanercept.json"],
        capture_output=True,
        text=True,
        timeout=30,
        cwd=repo_root(),
    )
    assert result.returncode == 0
    assert "COMPOSITE SCORE" in result.stdout


def test_cli_sensitivity_runs():
    result = subprocess.run(
        [sys.executable, "-m", "oglycan", "sensitivity", "sites/etanercept.json"],
        capture_output=True,
        text=True,
        timeout=30,
        cwd=repo_root(),
    )
    assert result.returncode == 0
    assert "Top sensitivities:" in result.stdout


def test_cli_compare_default_params():
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "oglycan",
            "compare",
            "sites/etanercept.json",
            "sites/epo.json",
        ],
        capture_output=True,
        text=True,
        timeout=30,
        cwd=repo_root(),
    )
    assert result.returncode == 0
    assert "etanercept" in result.stdout
    assert "epo" in result.stdout
    assert "0.9392" in result.stdout
    assert "0.9595" in result.stdout


def test_cli_compare_json():
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "oglycan",
            "compare",
            "sites/etanercept.json",
            "--json",
        ],
        capture_output=True,
        text=True,
        timeout=30,
        cwd=repo_root(),
    )
    assert result.returncode == 0
    payload = json.loads(result.stdout)
    assert isinstance(payload, list)
    assert len(payload) == 1


def test_cli_recommend_runs():
    result = subprocess.run(
        [sys.executable, "-m", "oglycan", "recommend", "sites/etanercept.json"],
        capture_output=True,
        text=True,
        timeout=30,
        cwd=repo_root(),
    )
    assert result.returncode == 0
    assert "COMPOSITE AFTER" in result.stdout


def test_cli_recommend_out(tmp_path):
    out_path = tmp_path / "recommend.json"
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "oglycan",
            "recommend",
            "sites/etanercept.json",
            "--out",
            str(out_path),
        ],
        capture_output=True,
        text=True,
        timeout=30,
        cwd=repo_root(),
    )
    assert result.returncode == 0
    assert out_path.exists()
    payload = json.loads(out_path.read_text())
    assert "recommended_params" in payload
    assert re.fullmatch(r"[0-9a-f]{64}", payload["source_catalog_sha256"])
    assert f"Recommended params written to {out_path}" in result.stdout


def test_cli_inspect_catalog():
    result = subprocess.run(
        [sys.executable, "-m", "oglycan", "inspect-catalog", "sites/etanercept.json"],
        capture_output=True,
        text=True,
        timeout=30,
        cwd=repo_root(),
    )
    assert result.returncode == 0
    assert "Etanercept" in result.stdout
    assert "13" in result.stdout
    assert "0.40" in result.stdout
    assert "0.90" in result.stdout


def test_cli_pipeline_writes_bundle(tmp_path):
    out_dir = tmp_path / "bundle"
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "oglycan",
            "pipeline",
            "sites/etanercept.json",
            "--out",
            str(out_dir),
        ],
        capture_output=True,
        text=True,
        timeout=30,
        cwd=repo_root(),
    )
    assert result.returncode == 0
    assert "Pipeline complete." in result.stdout

    for name in (
        "evaluate.json",
        "sensitivity.json",
        "recommend.json",
        "recommended_params.json",
        "method_card.md",
    ):
        path = out_dir / name
        assert path.exists()
        assert path.stat().st_size > 0


def test_cli_tune_warns_on_empty_catalog(tmp_path):
    catalog_path = tmp_path / "empty_catalog.json"
    catalog_path.write_text(
        json.dumps(
            {
                "glycoprotein": {"name": "EmptyCatalog"},
                "sites": [],
                "localization_threshold": 0.75,
            }
        ),
        encoding="utf-8",
    )

    result = subprocess.run(
        [sys.executable, "-m", "oglycan", "tune", str(catalog_path)],
        capture_output=True,
        text=True,
        timeout=30,
        cwd=repo_root(),
    )
    assert result.returncode == 0
    assert "⚠ WARNING: 0 scored sites in catalog." in result.stderr
    assert "Curate `sites[]` before relying on this score." in result.stderr


def test_cli_diff_catalog_detects_change(tmp_path):
    source_path = os.path.join(repo_root(), "sites", "etanercept.json")
    modified_path = tmp_path / "etanercept.modified.json"
    catalog = json.loads(open(source_path, encoding="utf-8").read())
    catalog["sites"][-1]["difficulty"] = 0.775
    modified_path.write_text(json.dumps(catalog, indent=2) + "\n", encoding="utf-8")

    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "oglycan",
            "diff-catalog",
            "sites/etanercept.json",
            str(modified_path),
        ],
        capture_output=True,
        text=True,
        timeout=30,
        cwd=repo_root(),
    )

    assert result.returncode == 0
    assert "delta" in result.stdout
    assert "T266" in result.stdout


def test_cli_diff_catalog_identical():
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "oglycan",
            "diff-catalog",
            "sites/etanercept.json",
            "sites/etanercept.json",
        ],
        capture_output=True,
        text=True,
        timeout=30,
        cwd=repo_root(),
    )

    assert result.returncode == 0
    assert result.stdout.strip() == "No changes."


def test_cli_suggest_catalog_runs(tmp_path):
    out_path = tmp_path / "toy.predicted.json"
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "oglycan",
            "suggest-catalog",
            "MSTPAAAPSTPAAAMSTPAA",
            "--name",
            "toy",
            "--out",
            str(out_path),
        ],
        capture_output=True,
        text=True,
        timeout=30,
        cwd=repo_root(),
    )
    assert result.returncode == 0
    assert out_path.exists()
    catalog = json.loads(out_path.read_text())
    assert catalog["glycoprotein"]["name"] == "toy"
    assert catalog["sites"] == []
    assert len(catalog["predicted_sites"]) == 6


def test_cli_ingest_pilot_o_pair_runs(tmp_path):
    pilot_path = tmp_path / "opair.tsv"
    pilot_path.write_text(
        "\n".join(
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
                        "Etanercept",
                        "AATTK",
                        "0.88",
                        "HexNAc(1)",
                        "[263-267]",
                        "4",
                    ]
                ),
            ]
        ),
        encoding="utf-8",
    )

    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "oglycan",
            "ingest-pilot",
            str(pilot_path),
            "--catalog",
            "sites/etanercept.json",
            "--format",
            "o_pair",
        ],
        capture_output=True,
        text=True,
        timeout=30,
        cwd=repo_root(),
    )

    assert result.returncode == 0
    assert "PILOT SOURCE: o_pair" in result.stdout
    assert "T266:" in result.stdout
