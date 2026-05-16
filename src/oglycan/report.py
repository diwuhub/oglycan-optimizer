"""Lean report generator: text + JSON output from evaluation results.

Replaces the 602-LOC report_source.py with ~70 LOC focused on
text summary and JSON export. HTML report generation is deferred
to v1.1 (not needed for MVP).
"""

import hashlib
import json
from datetime import datetime, timezone

from .method_card import _read_package_version


def text_report(result: dict) -> str:
    """Generate a text report from evaluation result dict."""
    lines = []
    w = 60

    lines.append("=" * w)
    lines.append("O-Glycan Site Localization Evaluator")
    lines.append("=" * w)
    lines.append(f"Reference: {result['molecule']} ({result['total_sites']} O-glycosites)")
    lines.append("")

    lines.append("Sub-model scores:")
    for name, val in result["sub_model_scores"].items():
        nice = name.replace("_", " ")
        lines.append(f"  {nice:<22} {val:.4f}")
    lines.append("")

    lines.append("Per-site localization:")
    for sr in result["site_results"]:
        status = "PASS" if sr["pass"] else "FAIL"
        forms = ", ".join(sr["core_types"])
        lines.append(
            f"  {sr['amino_acid']}{sr['position']}: {sr['confidence']:.4f} "
            f"[{status}]  (d={sr['difficulty']}, {forms})"
        )
    lines.append("")
    lines.append(f"Sites localized: {result['sites_localized']}/{result['total_sites']}")
    lines.append("")

    lines.append("Composite breakdown:")
    for name, val in result["composite_breakdown"].items():
        nice = name.replace("_", " ")
        lines.append(f"  {nice:<22} {val:.4f}")
    lines.append("")
    lines.append(f"COMPOSITE SCORE: {result['composite_score']:.4f}")

    return "\n".join(lines)


def _sha256_file(path: str) -> str:
    with open(path, "rb") as handle:
        return hashlib.sha256(handle.read()).hexdigest()


def json_report(
    result: dict,
    *,
    source_catalog: str | None = None,
    source_params: str | None = None,
    full_params: dict | None = None,
) -> str:
    """Generate a JSON report string from evaluation result dict."""
    output = {
        "generated_at": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "package_version": _read_package_version(),
        **result,
    }
    if source_catalog is not None:
        output["source_catalog_path"] = source_catalog
        output["source_catalog_sha256"] = _sha256_file(source_catalog)
    if source_params is not None:
        output["source_params_path"] = source_params
        output["source_params_sha256"] = _sha256_file(source_params)
    if full_params is not None:
        output["acquisition_params"] = full_params
    return json.dumps(output, indent=2)


def write_json(
    result: dict,
    path: str,
    *,
    source_catalog: str | None = None,
    source_params: str | None = None,
    full_params: dict | None = None,
) -> None:
    """Write JSON report to file."""
    with open(path, "w", encoding="utf-8") as f:
        f.write(
            json_report(
                result,
                source_catalog=source_catalog,
                source_params=source_params,
                full_params=full_params,
            )
        )
        f.write("\n")


def print_text(result: dict) -> None:
    """Print text report to stdout."""
    print(text_report(result))
