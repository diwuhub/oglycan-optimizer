"""CLI entry point for oglycan-optimizer.

Commands:
    serve           Launch the Streamlit UI
    compare         Compare evaluation results across one or more catalogs
    method-card     Write a markdown method card
    tune            Run evaluation and print text report
    uncertainty     Run evaluation and print composite score 95% CI
    sensitivity     Run signed local sensitivity analysis
    report          Run evaluation and write JSON report
    diff-catalog    Compare two site catalogs and print a diff
    inspect-catalog Print a one-screen catalog summary
    ingest-pilot    Load pilot evidence and print site summary
    recalibrate     Recalibrate site difficulties from pilot evidence
    recommend       Recommend next-run parameters under constraints
    pipeline        Run evaluate -> sensitivity -> recommend -> method-card
    suggest-catalog Build a starter catalog with predicted_sites candidates
    validate-sites  Validate a site catalog JSON

Usage:
    python -m oglycan serve
    python -m oglycan compare sites/etanercept.json sites/epo.json
    python -m oglycan method-card sites/etanercept.json --out method_card.md
    python -m oglycan tune sites/etanercept.json
    python -m oglycan tune etanercept
    python -m oglycan tune sites/etanercept.json --out result.json
    python -m oglycan uncertainty sites/etanercept.json
    python -m oglycan sensitivity sites/etanercept.json
    python -m oglycan tune sites/etanercept.json --params my_params.json
    python -m oglycan report sites/etanercept.json --out result.json
    python -m oglycan diff-catalog sites/etanercept.json sites/etanercept.json
    python -m oglycan inspect-catalog sites/etanercept.json
    python -m oglycan ingest-pilot pilot.tsv --catalog sites/etanercept.json
    python -m oglycan ingest-pilot pilot.tsv --catalog sites/etanercept.json --format o_pair
    python -m oglycan recalibrate pilot.json --catalog sites/etanercept.json --out sites/etanercept.recalibrated.json
    python -m oglycan recommend sites/etanercept.json --max-runtime 120
    python -m oglycan pipeline sites/etanercept.json --out bundle
    python -m oglycan suggest-catalog "MSTPAAAPSTPAAA" --name toy_protein
    python -m oglycan validate-sites sites/etanercept.json

Bare catalog names like `etanercept` resolve to `sites/etanercept.json`
relative to the repository root.
"""

import argparse
import copy
from datetime import datetime, timezone
import hashlib
import json
import os
import statistics
import subprocess
import sys

from .core import evaluate, load_params, load_site_catalog
from .method_card import _read_package_version, build_method_card
from .pilot import (
    load_pilot_canonical,
    load_pilot_msfragger,
    load_pilot_o_pair,
    validate_pilot,
)
from .predict import parse_fasta, suggest_catalog, validate_predicted_sites
from .recalibrate import recalibrate_difficulty, write_catalog
from .recommend import recommend_next_method
from .report import json_report, print_text, write_json
from .sensitivity import sensitivity
from .uncertainty import composite_with_ci


def _find_repo_root() -> str:
    """Find the repository root (contains examples/ and sites/ dirs).

    Tries cwd first (works for `python -m oglycan` from repo root),
    then falls back to package-relative path (works for PYTHONPATH=src).
    """
    cwd = os.getcwd()
    if os.path.isdir(os.path.join(cwd, "examples")) and os.path.isdir(os.path.join(cwd, "sites")):
        return cwd
    pkg_dir = os.path.dirname(os.path.abspath(__file__))
    repo_root = os.path.dirname(os.path.dirname(pkg_dir))
    if os.path.isdir(os.path.join(repo_root, "examples")):
        return repo_root
    return cwd


def _find_default_params() -> str:
    """Find default acquisition params relative to package root."""
    return os.path.join(_find_repo_root(), "examples", "default_acquisition_params.json")


def _find_web_app() -> str:
    """Find the Streamlit app relative to the package root."""
    pkg_dir = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(pkg_dir, "web", "app.py")


def _resolve_params_path(params_path: str | None) -> str:
    resolved_params_path = os.path.abspath(params_path or _find_default_params())
    if not os.path.exists(resolved_params_path):
        print(f"ERROR: params file not found: {resolved_params_path}", file=sys.stderr)
        sys.exit(1)
    return resolved_params_path


def _resolve_catalog_path(value: str) -> str:
    """
    If `value` is a bare name (no '/' or '\\' and doesn't end in .json),
    try resolving as sites/<value>.json relative to the repo root.
    Otherwise return verbatim.
    """
    if "/" in value or "\\" in value or value.lower().endswith(".json"):
        return value
    return os.path.join(_find_repo_root(), "sites", f"{value}.json")


def _load_catalog_and_params(
    sites_path: str,
    params_path: str | None,
) -> tuple[dict, dict, str, str]:
    resolved_sites_path = _resolve_catalog_path(sites_path)
    catalog = load_site_catalog(resolved_sites_path)
    resolved_params_path = _resolve_params_path(params_path)
    params = load_params(resolved_params_path)
    return catalog, params, resolved_sites_path, resolved_params_path


def _warn_if_no_scored_sites(catalog: dict, *, catalog_label: str | None = None) -> None:
    if not catalog.get("sites"):
        lines = ["⚠ WARNING: 0 scored sites in catalog."]
        if catalog_label:
            lines.append(f"  Catalog: {catalog_label}")
        lines.extend(
            [
                "  Composite reflects method-only dimensions (spectral quality, bio plausibility, glycan diversity). Localization and",
                "  coverage components are zero. Curate `sites[]` before relying on this score.",
            ]
        )
        print("\n".join(lines), file=sys.stderr)


def _utc_timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _sha256_file(path: str) -> str:
    with open(path, "rb") as handle:
        return hashlib.sha256(handle.read()).hexdigest()


def _write_json_payload(path: str, payload: dict) -> None:
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2)
        handle.write("\n")


def _catalog_name(catalog: dict) -> str:
    return catalog.get("glycoprotein", {}).get("name", "unknown")


def _catalog_label(catalog_path: str) -> str:
    return os.path.splitext(os.path.basename(catalog_path))[0]


def _ensure_pilot_matches_catalog(pilot: dict, catalog: dict) -> None:
    pilot_name = pilot.get("glycoprotein_name", "")
    catalog_name = _catalog_name(catalog)
    if pilot_name.lower() != catalog_name.lower():
        print(
            f"ERROR: pilot glycoprotein {pilot_name!r} does not match catalog {catalog_name!r}",
            file=sys.stderr,
        )
        sys.exit(1)


def _resolve_pilot_format(path: str, requested_format: str | None) -> str:
    if requested_format:
        return requested_format
    ext = os.path.splitext(path)[1].lower()
    if ext == ".json":
        return "canonical"
    if ext == ".tsv":
        return "msfragger"
    print(
        "ERROR: could not infer pilot format from extension; use --format canonical, --format msfragger, or --format o_pair",
        file=sys.stderr,
    )
    sys.exit(1)


def _load_pilot(path: str, pilot_format: str, catalog: dict) -> dict:
    if pilot_format == "canonical":
        pilot = load_pilot_canonical(path)
    elif pilot_format == "msfragger":
        pilot = load_pilot_msfragger(path, glycoprotein_name=_catalog_name(catalog))
    elif pilot_format == "o_pair":
        pilot = load_pilot_o_pair(path, glycoprotein_name=_catalog_name(catalog))
    else:
        print(f"ERROR: unsupported pilot format: {pilot_format}", file=sys.stderr)
        sys.exit(1)

    errors = validate_pilot(pilot)
    if errors:
        print("ERROR: invalid pilot result:", file=sys.stderr)
        for error in errors:
            print(f"  - {error}", file=sys.stderr)
        sys.exit(1)

    _ensure_pilot_matches_catalog(pilot, catalog)
    return pilot


def _print_pilot_summary(pilot: dict) -> None:
    print(f"PILOT SOURCE: {pilot['source']}")
    print(f"GLYCOPROTEIN: {pilot['glycoprotein_name']}")
    print(f"SITES: {len(pilot['sites'])}")
    for site in sorted(pilot["sites"], key=lambda item: (item["pos"], item["aa"])):
        loc = (
            "None"
            if site["observed_localization"] is None
            else f"{site['observed_localization']:.4f}"
        )
        glycoforms = ", ".join(site["observed_glycoforms"]) or "(none)"
        print(
            f"{site['aa']}{site['pos']}: "
            f"localization={loc} n_spectra={site['n_spectra']} glycoforms={glycoforms}"
        )


def _load_sequence_input(sequence_or_path: str) -> str:
    ext = os.path.splitext(sequence_or_path)[1].lower()
    if os.path.isfile(sequence_or_path) and ext in {".fa", ".fasta", ".txt"}:
        with open(sequence_or_path) as handle:
            records = parse_fasta(handle.read())
        if not records:
            print(
                f"ERROR: no FASTA records found in {sequence_or_path}",
                file=sys.stderr,
            )
            sys.exit(1)
        return records[0][1]
    return sequence_or_path


def _catalog_validation_errors(catalog: dict) -> list[str]:
    sites = catalog.get("sites", [])
    errors = []
    if not isinstance(sites, list):
        return ["'sites' must be a list"]

    for i, site in enumerate(sites):
        if not isinstance(site, dict):
            errors.append(f"site[{i}]: must be an object")
            continue
        for field in ("pos", "aa", "core_types", "difficulty"):
            if field not in site:
                errors.append(f"site[{i}]: missing '{field}'")
        if "pos" in site:
            pos = site["pos"]
            if isinstance(pos, bool) or not isinstance(pos, int):
                errors.append(f"site[{i}]: pos must be an integer")
            elif pos < 1:
                errors.append(f"site[{i}]: pos must be >= 1")
        if "aa" in site:
            aa = str(site["aa"]).strip().upper()
            if aa not in {"S", "T"}:
                errors.append(f"site[{i}]: aa must be 'S' or 'T'")
        if "core_types" in site:
            core_types = site["core_types"]
            if not isinstance(core_types, list):
                errors.append(f"site[{i}]: core_types must be a list")
            elif not [item for item in core_types if str(item).strip()]:
                errors.append(f"site[{i}]: core_types must be non-empty")
        if "difficulty" in site:
            try:
                difficulty = float(site["difficulty"])
            except (TypeError, ValueError):
                errors.append(f"site[{i}]: difficulty must be numeric")
            else:
                if not (0.0 <= difficulty <= 1.0):
                    errors.append(f"site[{i}]: difficulty {difficulty} out of [0, 1]")

    predicted_sites = catalog.get("predicted_sites")
    if predicted_sites is not None:
        errors.extend(validate_predicted_sites(predicted_sites))
    return errors


def _write_evaluation_report_json(
    result: dict,
    *,
    out_path: str,
    source_catalog: str,
    source_params: str,
    full_params: dict,
) -> None:
    write_json(
        result,
        out_path,
        source_catalog=source_catalog,
        source_params=source_params,
        full_params=full_params,
    )


def _compute_catalog_diff(old_catalog: dict, new_catalog: dict) -> dict:
    old_sites = old_catalog.get("sites") or []
    new_sites = new_catalog.get("sites") or []
    old_by_pos = {site["pos"]: site for site in old_sites}
    new_by_pos = {site["pos"]: site for site in new_sites}

    sites_changed = []
    for old_site in old_sites:
        pos = old_site["pos"]
        new_site = new_by_pos.get(pos)
        if new_site is None:
            continue
        difficulty_old = float(old_site["difficulty"])
        difficulty_new = float(new_site["difficulty"])
        core_types_old = old_site.get("core_types", [])
        core_types_new = new_site.get("core_types", [])
        if difficulty_old != difficulty_new or sorted(core_types_old) != sorted(core_types_new):
            sites_changed.append(
                {
                    "pos": pos,
                    "aa": new_site.get("aa", old_site.get("aa", "")),
                    "difficulty_old": difficulty_old,
                    "difficulty_new": difficulty_new,
                    "core_types_old": list(core_types_old),
                    "core_types_new": list(core_types_new),
                }
            )

    sites_added = [copy.deepcopy(site) for site in new_sites if site["pos"] not in old_by_pos]
    sites_removed = [copy.deepcopy(site) for site in old_sites if site["pos"] not in new_by_pos]

    return {
        "catalog_name_old": _catalog_name(old_catalog),
        "catalog_name_new": _catalog_name(new_catalog),
        "sites_changed": sites_changed,
        "sites_added": sites_added,
        "sites_removed": sites_removed,
        "threshold_old": old_catalog.get("localization_threshold"),
        "threshold_new": new_catalog.get("localization_threshold"),
    }


def _has_catalog_diff(diff: dict) -> bool:
    return bool(
        diff["sites_changed"]
        or diff["sites_added"]
        or diff["sites_removed"]
        or diff["threshold_old"] != diff["threshold_new"]
    )


def _format_catalog_number(value: float | int | None) -> str:
    if value is None:
        return "None"
    formatted = f"{float(value):.3f}"
    if formatted.endswith("0"):
        formatted = formatted[:-1]
    return formatted


def _format_core_types(core_types: list[str]) -> str:
    return f"[{', '.join(core_types)}]"


def _format_site_summary(site: dict) -> str:
    parts = [
        f"{site.get('aa', '?')}{site.get('pos', '?')}",
        f"difficulty={_format_catalog_number(site['difficulty'])}",
        f"core_types={_format_core_types(site.get('core_types', []))}",
    ]
    if site.get("source"):
        parts.append(f"source={site['source']}")
    return ", ".join(parts)


def _format_catalog_diff_text(old_catalog: dict, new_catalog: dict, diff: dict) -> str:
    lines = [
        f"Catalog diff: {diff['catalog_name_old']} -> {diff['catalog_name_new']}",
        "",
        "Sites changed:",
    ]
    if diff["sites_changed"]:
        for item in diff["sites_changed"]:
            label = f"{item['aa']}{item['pos']}"
            if item["difficulty_old"] == item["difficulty_new"]:
                lines.append(f"  {label}  difficulty: unchanged")
            else:
                delta = item["difficulty_new"] - item["difficulty_old"]
                lines.append(
                    f"  {label}  difficulty: {_format_catalog_number(item['difficulty_old'])} -> "
                    f"{_format_catalog_number(item['difficulty_new'])} "
                    f"(delta {_format_catalog_number(delta)})"
                )
            if sorted(item["core_types_old"]) == sorted(item["core_types_new"]):
                lines.append("        core_types: unchanged")
            else:
                lines.append(
                    "        core_types: changed: "
                    f"{_format_core_types(item['core_types_old'])} -> "
                    f"{_format_core_types(item['core_types_new'])}"
                )
    else:
        lines.append("  (none)")

    lines.extend(["", "Sites added:"])
    if diff["sites_added"]:
        for site in diff["sites_added"]:
            lines.append(f"  {_format_site_summary(site)}")
    else:
        lines.append("  (none)")

    lines.extend(["", "Sites removed:"])
    if diff["sites_removed"]:
        for site in diff["sites_removed"]:
            lines.append(f"  {_format_site_summary(site)}")
    else:
        lines.append("  (none)")

    old_predicted = len(old_catalog.get("predicted_sites") or [])
    new_predicted = len(new_catalog.get("predicted_sites") or [])
    predicted_status = "(unchanged)" if old_predicted == new_predicted else "(changed)"
    threshold_status = (
        "(unchanged)"
        if diff["threshold_old"] == diff["threshold_new"]
        else "(changed)"
    )
    lines.extend(
        [
            "",
            f"Predicted sites: {old_predicted} -> {new_predicted} {predicted_status}",
            "Localization threshold: "
            f"{_format_catalog_number(diff['threshold_old'])} -> "
            f"{_format_catalog_number(diff['threshold_new'])} {threshold_status}",
        ]
    )
    return "\n".join(lines)


def _recommend_constraints(args) -> dict:
    constraints = {}
    if args.max_runtime is not None:
        constraints["max_runtime_min"] = args.max_runtime
    if args.min_composite is not None:
        constraints["min_composite"] = args.min_composite
    if args.instrument_class is not None:
        constraints["instrument_class"] = args.instrument_class
    return constraints


def _recommend_output(
    result: dict,
    *,
    source_catalog: str,
    source_params: str,
    constraints: dict,
) -> dict:
    return {
        "generated_at": _utc_timestamp(),
        "package_version": _read_package_version(),
        "source_catalog": source_catalog,
        "source_catalog_sha256": _sha256_file(source_catalog),
        "source_params": source_params,
        "source_params_sha256": _sha256_file(source_params),
        "constraints": constraints,
        "composite_before": result["composite_before"],
        "composite_after": result["composite_after"],
        "improvement": result["improvement"],
        "changes": result["changes"],
        "constraint_violations": result["constraint_violations"],
        "recommended_params": result["recommended_params"],
    }


def _sensitivity_summary_delta(row: dict) -> float:
    if abs(row["delta_down"]) > abs(row["delta_up"]):
        return row["delta_down"]
    return row["delta_up"]


def _compare_row(
    catalog_path: str,
    catalog: dict,
    result: dict,
    *,
    top_sensitivity_rows: list[dict] | None = None,
) -> dict:
    row = {
        "catalog": _catalog_label(catalog_path),
        "catalog_path": catalog_path,
        "glycoprotein_name": _catalog_name(catalog),
        "sites": result["total_sites"],
        "localized_sites": result["sites_localized"],
        "localized_label": f"{result['sites_localized']}/{result['total_sites']}",
        "composite": float(result["composite_score"]),
        "localization_confidence": float(
            result["composite_breakdown"]["localization_confidence"]
        ),
        "sequence_coverage": float(result["composite_breakdown"]["sequence_coverage"]),
        "spectral_quality": float(result["composite_breakdown"]["spectral_quality"]),
    }
    if top_sensitivity_rows is not None:
        row["top_sensitivity_movers"] = [
            {
                "param_path": item["param_path"],
                "label": item["param_path"].split(".")[-1],
                "delta": _sensitivity_summary_delta(item),
            }
            for item in top_sensitivity_rows
        ]
    return row


def _format_compare_text(
    rows: list[dict],
    *,
    params_label: str,
    include_sensitivity: bool,
) -> str:
    catalog_width = max(len("Catalog"), *(len(row["catalog"]) for row in rows))
    sites_width = max(len("Sites"), *(len(str(row["sites"])) for row in rows))
    localized_width = max(
        len("Localized"),
        *(len(row["localized_label"]) for row in rows),
    )

    lines = [
        f"Cross-protein comparison ({len(rows)} catalogs, {params_label})",
        "",
        f"{'Catalog':<{catalog_width}}  "
        f"{'Sites':>{sites_width}}  "
        f"{'Localized':>{localized_width}}  "
        f"{'Composite':>9}  "
        f"{'LocConf':>7}  "
        f"{'SeqCov':>6}  "
        f"{'SpecQual':>8}",
    ]

    for row in rows:
        lines.append(
            f"{row['catalog']:<{catalog_width}}  "
            f"{row['sites']:>{sites_width}}  "
            f"{row['localized_label']:>{localized_width}}  "
            f"{row['composite']:>9.4f}  "
            f"{row['localization_confidence']:>7.4f}  "
            f"{row['sequence_coverage']:>6.2f}  "
            f"{row['spectral_quality']:>8.4f}"
        )

    if include_sensitivity:
        lines.extend(["", "Top sensitivity movers per catalog:"])
        for row in rows:
            movers = row.get("top_sensitivity_movers") or []
            summary = (
                ", ".join(
                    f"{item['label']} ({item['delta']:.4f})" for item in movers
                )
                if movers
                else "none"
            )
            lines.append(f"  {row['catalog']:<{catalog_width}}  {summary}")

    return "\n".join(lines)


def cmd_tune(args):
    catalog, params, resolved_sites_path, resolved_params_path = _load_catalog_and_params(
        args.sites,
        args.params,
    )
    _warn_if_no_scored_sites(catalog)
    result = evaluate(catalog, params)
    if args.out:
        _write_evaluation_report_json(
            result,
            out_path=args.out,
            source_catalog=resolved_sites_path,
            source_params=resolved_params_path,
            full_params=params,
        )
    if args.quiet:
        return
    if args.json:
        print(
            json_report(
                result,
                source_catalog=resolved_sites_path,
                source_params=resolved_params_path,
                full_params=params,
            )
        )
    else:
        print_text(result)


def cmd_serve(_args):
    try:
        import streamlit  # noqa: F401
    except ImportError:
        print(
            "ERROR: streamlit is not installed. Install optional web deps with `pip install .[web]`.",
            file=sys.stderr,
        )
        sys.exit(1)

    try:
        subprocess.run(["streamlit", "run", _find_web_app()], check=False)
    except FileNotFoundError:
        print(
            "ERROR: `streamlit` executable not found. Install with `pip install .[web]`.",
            file=sys.stderr,
        )
        sys.exit(1)


def cmd_method_card(args):
    catalog, params, _resolved_sites_path, _resolved_params_path = _load_catalog_and_params(
        args.sites,
        args.params,
    )
    _warn_if_no_scored_sites(catalog)
    result = evaluate(catalog, params)
    out_path = args.out or "method_card.md"
    with open(out_path, "w", encoding="utf-8") as handle:
        handle.write(build_method_card(catalog, params, result))
    print(f"Method card written to {out_path}")


def cmd_uncertainty(args):
    catalog, params, _resolved_sites_path, _resolved_params_path = _load_catalog_and_params(
        args.sites,
        args.params,
    )
    _warn_if_no_scored_sites(catalog)
    result = composite_with_ci(catalog, params)
    print(
        "COMPOSITE SCORE: "
        f"{result['point']:.4f} "
        f"(95% CI {result['ci_low']:.4f}-{result['ci_high']:.4f}, "
        f"n={result['n_samples']})"
    )


def cmd_sensitivity(args):
    catalog, params, _resolved_sites_path, _resolved_params_path = _load_catalog_and_params(
        args.sites,
        args.params,
    )
    _warn_if_no_scored_sites(catalog)
    rows = sensitivity(catalog, params)
    if args.json:
        print(json.dumps(rows, indent=2))
        return
    print("Top sensitivities:")
    print(f"{'param_path':<42} {'kind':<8} {'current':>12} {'delta_up':>10} {'delta_down':>10} {'abs':>10}")
    for row in rows[:10]:
        print(
            f"{row['param_path']:<42} "
            f"{row['kind']:<8} "
            f"{str(row['current_value']):>12} "
            f"{row['delta_up']:>10.4f} "
            f"{row['delta_down']:>10.4f} "
            f"{row['abs_effect']:>10.4f}"
        )


def cmd_report(args):
    catalog, params, resolved_sites_path, resolved_params_path = _load_catalog_and_params(
        args.sites,
        args.params,
    )
    _warn_if_no_scored_sites(catalog)
    result = evaluate(catalog, params)
    out_path = args.out or "report.json"
    _write_evaluation_report_json(
        result,
        out_path=out_path,
        source_catalog=resolved_sites_path,
        source_params=resolved_params_path,
        full_params=params,
    )
    print(f"JSON report written to {out_path}")


def cmd_diff_catalog(args):
    try:
        old_catalog = load_site_catalog(_resolve_catalog_path(args.old_catalog))
        new_catalog = load_site_catalog(_resolve_catalog_path(args.new_catalog))
    except (json.JSONDecodeError, FileNotFoundError, OSError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)

    diff = _compute_catalog_diff(old_catalog, new_catalog)
    if args.json:
        print(json.dumps(diff, indent=2))
    elif _has_catalog_diff(diff):
        print(_format_catalog_diff_text(old_catalog, new_catalog, diff))
    else:
        print("No changes.")


def cmd_inspect_catalog(args):
    try:
        catalog = load_site_catalog(_resolve_catalog_path(args.sites))
    except (json.JSONDecodeError, FileNotFoundError, OSError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)

    errors = _catalog_validation_errors(catalog)
    if errors:
        print(f"ERROR: invalid catalog: {errors[0]}", file=sys.stderr)
        sys.exit(1)

    glycoprotein = catalog.get("glycoprotein", {})
    sites = catalog.get("sites", [])
    predicted_sites = catalog.get("predicted_sites") or []
    n_glycosites = catalog.get("n_glycosites") or []
    threshold = catalog.get("localization_threshold")

    def display(value) -> str:
        return "(not set)" if value in (None, "") else str(value)

    if sites:
        difficulties = [float(site["difficulty"]) for site in sites]
        scored_sites_summary = (
            f"{len(sites)}  (difficulty min={min(difficulties):.2f}, "
            f"median={statistics.median(difficulties):.2f}, max={max(difficulties):.2f})"
        )
    else:
        scored_sites_summary = "0  (difficulty min=n/a, median=n/a, max=n/a)"

    immunogenicity_counts = {}
    for site in sites:
        for flag in site.get("immunogenicity_flags", []):
            immunogenicity_counts[flag] = immunogenicity_counts.get(flag, 0) + 1
    if immunogenicity_counts:
        immunogenicity_summary = ", ".join(
            f"{flag} ({count} sites)"
            for flag, count in sorted(
                immunogenicity_counts.items(),
                key=lambda item: (-item[1], item[0]),
            )
        )
    else:
        immunogenicity_summary = "none"

    threshold_value = display(threshold)
    if isinstance(threshold, (int, float)):
        threshold_value = f"{threshold:.2f}"

    print(f"Catalog: {display(glycoprotein.get('name'))}")
    print(f"  {'UniProt:':<19} {display(glycoprotein.get('uniprot'))}")
    print(f"  {'Reference:':<19} {display(glycoprotein.get('reference'))}")
    print(f"  {'Scored sites:':<19} {scored_sites_summary}")
    print(f"  {'Predicted sites:':<19} {len(predicted_sites)}")
    print(f"  {'N-glycosites:':<19} {len(n_glycosites)}")
    print(f"  {'Localization threshold:':<19} {threshold_value}")
    print(f"  {'Immunogenicity flags present:':<19} {immunogenicity_summary}")


def cmd_ingest_pilot(args):
    catalog = load_site_catalog(_resolve_catalog_path(args.catalog))
    pilot_format = _resolve_pilot_format(args.pilot, args.format)
    pilot = _load_pilot(args.pilot, pilot_format, catalog)
    _print_pilot_summary(pilot)


def cmd_recalibrate(args):
    catalog = load_site_catalog(_resolve_catalog_path(args.catalog))
    pilot = _load_pilot(args.pilot, "canonical", catalog)
    recalibrated = recalibrate_difficulty(catalog, pilot)
    write_catalog(recalibrated, args.out)
    meta = recalibrated["recalibration_metadata"]
    print(f"Recalibrated catalog written to {args.out}")
    print(
        "Updated={n_sites_updated} Unobserved={n_sites_unobserved} Added={n_sites_added}".format(
            **meta
        )
    )


def cmd_recommend(args):
    catalog, params, resolved_sites_path, resolved_params_path = _load_catalog_and_params(
        args.sites,
        args.params,
    )
    _warn_if_no_scored_sites(catalog)
    constraints = _recommend_constraints(args)

    result = recommend_next_method(
        catalog,
        params,
        constraints if constraints else None,
    )
    payload = _recommend_output(
        result,
        source_catalog=resolved_sites_path,
        source_params=resolved_params_path,
        constraints=constraints,
    )
    if args.out:
        _write_json_payload(args.out, payload)
    if args.json:
        print(json.dumps(payload, indent=2))
        return
    print(f"COMPOSITE BEFORE: {result['composite_before']:.4f}")
    print(f"COMPOSITE AFTER:  {result['composite_after']:.4f}")
    print(f"COMPOSITE Δ:      {result['improvement']:.4f}")
    print("CHANGES:")
    if result["changes"]:
        for change in result["changes"]:
            print(
                f"  - {change['path']}: {change['before']} -> {change['after']}"
            )
    else:
        print("  - none")
    print("CONSTRAINT VIOLATIONS:")
    if result["constraint_violations"]:
        for violation in result["constraint_violations"]:
            print(f"  - {violation}")
    else:
        print("  - none")
    if args.out:
        print(f"Recommended params written to {args.out}")


def cmd_pipeline(args):
    catalog, params, resolved_sites_path, resolved_params_path = _load_catalog_and_params(
        args.sites,
        args.params,
    )
    _warn_if_no_scored_sites(catalog)
    constraints = _recommend_constraints(args)

    out_dir = os.path.abspath(args.out)
    os.makedirs(out_dir, exist_ok=True)

    evaluation = evaluate(catalog, params)
    evaluate_path = os.path.join(out_dir, "evaluate.json")
    write_json(
        evaluation,
        evaluate_path,
        source_catalog=resolved_sites_path,
        source_params=resolved_params_path,
        full_params=params,
    )

    sensitivity_rows = sensitivity(catalog, params)
    sensitivity_path = os.path.join(out_dir, "sensitivity.json")
    _write_json_payload(
        sensitivity_path,
        {
            "generated_at": _utc_timestamp(),
            "package_version": _read_package_version(),
            "sensitivity": sensitivity_rows[:10],
            "full_ranking": sensitivity_rows,
        },
    )

    recommendation = recommend_next_method(
        catalog,
        params,
        constraints if constraints else None,
    )
    recommend_path = os.path.join(out_dir, "recommend.json")
    _write_json_payload(
        recommend_path,
        _recommend_output(
            recommendation,
            source_catalog=resolved_sites_path,
            source_params=resolved_params_path,
            constraints=constraints,
        ),
    )

    recommended_params = recommendation["recommended_params"]
    recommended_params_path = os.path.join(out_dir, "recommended_params.json")
    _write_json_payload(recommended_params_path, recommended_params)

    recommended_evaluation = evaluate(catalog, recommended_params)
    method_card = build_method_card(catalog, recommended_params, recommended_evaluation)
    method_card_path = os.path.join(out_dir, "method_card.md")
    with open(method_card_path, "w", encoding="utf-8") as handle:
        handle.write(method_card)

    top_mover = "none"
    if sensitivity_rows:
        top_row = sensitivity_rows[0]
        top_mover = f"{top_row['param_path']} ({_sensitivity_summary_delta(top_row):.4f})"

    print(f"Pipeline complete. Wrote to {out_dir}/")
    print(
        f"  evaluate.json            composite {evaluation['composite_score']:.4f} (with defaults)"
    )
    print(f"  sensitivity.json         top mover: {top_mover}")
    print(
        "  recommend.json           "
        f"composite {recommendation['composite_after']:.4f} after constraints "
        f"(delta {recommendation['improvement']:.4f})"
    )
    print("  recommended_params.json  full params dict for downstream re-use")
    print(
        f"  method_card.md           {len(method_card.splitlines())} lines, ready for manual transcription"
    )


def cmd_compare(args):
    resolved_params_path = _resolve_params_path(args.params)
    params = load_params(resolved_params_path)
    rows = []
    warned_catalogs = set()

    for catalog_arg in args.catalogs:
        resolved_catalog_path = _resolve_catalog_path(catalog_arg)
        catalog = load_site_catalog(resolved_catalog_path)
        catalog_label = _catalog_label(resolved_catalog_path)
        if not catalog.get("sites") and catalog_label not in warned_catalogs:
            _warn_if_no_scored_sites(catalog, catalog_label=catalog_label)
            warned_catalogs.add(catalog_label)
        result = evaluate(catalog, params)
        top_sensitivity_rows = None
        if args.sensitivity:
            top_sensitivity_rows = sensitivity(catalog, params)[:3]
        rows.append(
            _compare_row(
                resolved_catalog_path,
                catalog,
                result,
                top_sensitivity_rows=top_sensitivity_rows,
            )
        )

    if args.json:
        print(json.dumps(rows, indent=2))
        return

    params_label = "default params" if args.params is None else f"params: {resolved_params_path}"
    print(
        _format_compare_text(
            rows,
            params_label=params_label,
            include_sensitivity=args.sensitivity,
        )
    )


def cmd_suggest_catalog(args):
    sequence = _load_sequence_input(args.sequence_or_fasta_path)
    out_path = args.out or f"{args.name}.predicted.json"
    try:
        catalog = suggest_catalog(
            sequence=sequence,
            protein_name=args.name,
            predictor=args.predictor,
        )
    except (NotImplementedError, RuntimeError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)

    with open(out_path, "w") as handle:
        json.dump(catalog, handle, indent=2)
        handle.write("\n")
    print(f"Starter catalog written to {out_path}")


def cmd_validate(args):
    try:
        catalog = load_site_catalog(_resolve_catalog_path(args.sites))
    except (json.JSONDecodeError, FileNotFoundError) as e:
        print(f"FAIL: {e}", file=sys.stderr)
        sys.exit(1)

    sites = catalog.get("sites", [])
    predicted_sites = catalog.get("predicted_sites")
    errors = _catalog_validation_errors(catalog)

    if errors:
        print(f"FAIL: {len(errors)} error(s):")
        for e in errors:
            print(f"  - {e}")
        sys.exit(1)
    else:
        name = catalog.get("glycoprotein", {}).get("name", "unknown")
        n_predicted = len(predicted_sites) if predicted_sites is not None else 0
        if predicted_sites is None:
            print(f"OK: {name}, {len(sites)} sites, all valid.")
        else:
            print(
                f"OK: {name}, {len(sites)} sites, {n_predicted} predicted_sites, all valid."
            )


def main():
    parser = argparse.ArgumentParser(
        prog="oglycan",
        description="O-glycan site localization optimizer for O-glycopeptide MS"
    )
    sub = parser.add_subparsers(dest="command")

    sub.add_parser("serve", help="Launch the Streamlit UI")

    p_compare = sub.add_parser("compare", help="Compare evaluation results across one or more catalogs")
    p_compare.add_argument("catalogs", nargs="+", help="Catalog path(s) or bare catalog names")
    p_compare.add_argument("--params", help="Path to acquisition params JSON")
    p_compare.add_argument("--json", action="store_true", help="Emit structured JSON rows")
    p_compare.add_argument(
        "--sensitivity",
        action="store_true",
        help="Include top-3 sensitivity movers per catalog",
    )

    p_method_card = sub.add_parser("method-card", help="Write a markdown method card")
    p_method_card.add_argument("sites", help="Path to site catalog JSON")
    p_method_card.add_argument("--params", help="Path to acquisition params JSON")
    p_method_card.add_argument("--out", help="Output markdown path (default: method_card.md)")

    p_tune = sub.add_parser("tune", help="Run evaluation and print text report")
    p_tune.add_argument("sites", help="Path to site catalog JSON")
    p_tune.add_argument("--params", help="Path to acquisition params JSON (default: examples/default_acquisition_params.json)")
    p_tune.add_argument("--out", help="Write the JSON report to this path in addition to the text report")
    p_tune.add_argument("--json", action="store_true", help="Emit the JSON report to stdout instead of text")
    p_tune.add_argument("--quiet", action="store_true", help="Suppress the text report output; useful with --out")

    p_uncertainty = sub.add_parser("uncertainty", help="Run evaluation and print composite score 95 percent CI")
    p_uncertainty.add_argument("sites", help="Path to site catalog JSON")
    p_uncertainty.add_argument("--params", help="Path to acquisition params JSON")

    p_sensitivity = sub.add_parser("sensitivity", help="Run signed local sensitivity analysis")
    p_sensitivity.add_argument("sites", help="Path to site catalog JSON")
    p_sensitivity.add_argument("--params", help="Path to acquisition params JSON")
    p_sensitivity.add_argument("--json", action="store_true", help="Emit the full ranked sensitivity list as JSON")

    p_report = sub.add_parser(
        "report",
        help="Run evaluation and write JSON report (prefer 'tune --out'; retained for backward compatibility)",
    )
    p_report.add_argument("sites", help="Path to site catalog JSON")
    p_report.add_argument("--params", help="Path to acquisition params JSON")
    p_report.add_argument("--out", help="Output JSON path (default: report.json)")

    p_diff = sub.add_parser("diff-catalog", help="Compare two site catalogs and print a diff")
    p_diff.add_argument("old_catalog", help="Path to the original site catalog JSON")
    p_diff.add_argument("new_catalog", help="Path to the updated site catalog JSON")
    p_diff.add_argument("--json", action="store_true", help="Emit the structured diff as JSON")

    p_inspect = sub.add_parser("inspect-catalog", help="Print a one-screen catalog summary")
    p_inspect.add_argument("sites", help="Path to site catalog JSON")

    p_ingest = sub.add_parser("ingest-pilot", help="Load pilot evidence and print a site summary")
    p_ingest.add_argument(
        "pilot",
        help="Path to canonical pilot JSON or pilot TSV (MSFragger/O-Pair)",
    )
    p_ingest.add_argument("--catalog", required=True, help="Path to site catalog JSON")
    p_ingest.add_argument(
        "--format",
        choices=["canonical", "msfragger", "o_pair"],
        help=(
            "Pilot input format; defaults to extension-based detection "
            "(.json=canonical, .tsv=msfragger). Use --format o_pair for O-Pair TSV."
        ),
    )

    p_recal = sub.add_parser("recalibrate", help="Write a recalibrated site catalog from canonical pilot JSON")
    p_recal.add_argument("pilot", help="Path to canonical pilot JSON")
    p_recal.add_argument("--catalog", required=True, help="Path to site catalog JSON")
    p_recal.add_argument("--out", required=True, help="Output path for recalibrated catalog JSON")

    p_recommend = sub.add_parser("recommend", help="Recommend next-run parameters under optional constraints")
    p_recommend.add_argument("sites", help="Path to site catalog JSON")
    p_recommend.add_argument("--params", help="Path to acquisition params JSON")
    p_recommend.add_argument("--max-runtime", type=float, help="Maximum allowed runtime in minutes")
    p_recommend.add_argument("--min-composite", type=float, help="Target composite score floor")
    p_recommend.add_argument("--out", help="Write recommended params bundle JSON to this path")
    p_recommend.add_argument("--json", action="store_true", help="Emit the full recommendation JSON to stdout")
    p_recommend.add_argument(
        "--instrument-class",
        choices=["orbitrap", "tims"],
        help="Instrument family to target",
    )

    p_pipeline = sub.add_parser("pipeline", help="Run evaluate, sensitivity, recommend, and method-card in one bundle")
    p_pipeline.add_argument("sites", help="Path to site catalog JSON")
    p_pipeline.add_argument("--params", help="Path to acquisition params JSON")
    p_pipeline.add_argument("--max-runtime", type=float, help="Maximum allowed runtime in minutes")
    p_pipeline.add_argument("--min-composite", type=float, help="Target composite score floor")
    p_pipeline.add_argument(
        "--instrument-class",
        choices=["orbitrap", "tims"],
        help="Instrument family to target",
    )
    p_pipeline.add_argument("--out", required=True, help="Bundle output directory")

    p_suggest = sub.add_parser(
        "suggest-catalog",
        help="Build a starter catalog with predicted_sites candidates",
        epilog=(
            "--predictor netoglyc makes a live HTTPS call to "
            "services.healthtech.dtu.dk. The service is intermittent; on failure, "
            "re-run later or use --predictor scan_st."
        ),
    )
    p_suggest.add_argument(
        "sequence_or_fasta_path",
        help="Raw protein sequence or FASTA path (.fa/.fasta/.txt)",
    )
    p_suggest.add_argument(
        "--name",
        required=True,
        help="Protein name for the starter catalog",
    )
    p_suggest.add_argument(
        "--predictor",
        choices=["scan_st", "netoglyc", "stackoglypred"],
        default="scan_st",
        help="Predictor to use; scan_st emits every S/T as a candidate for review",
    )
    p_suggest.add_argument(
        "--out",
        help="Output JSON path (default: <name>.predicted.json)",
    )

    p_val = sub.add_parser("validate-sites", help="Validate a site catalog JSON")
    p_val.add_argument("sites", help="Path to site catalog JSON")

    args = parser.parse_args()
    if args.command == "serve":
        cmd_serve(args)
    elif args.command == "compare":
        cmd_compare(args)
    elif args.command == "method-card":
        cmd_method_card(args)
    elif args.command == "tune":
        cmd_tune(args)
    elif args.command == "uncertainty":
        cmd_uncertainty(args)
    elif args.command == "sensitivity":
        cmd_sensitivity(args)
    elif args.command == "report":
        cmd_report(args)
    elif args.command == "diff-catalog":
        cmd_diff_catalog(args)
    elif args.command == "inspect-catalog":
        cmd_inspect_catalog(args)
    elif args.command == "ingest-pilot":
        cmd_ingest_pilot(args)
    elif args.command == "recalibrate":
        cmd_recalibrate(args)
    elif args.command == "recommend":
        cmd_recommend(args)
    elif args.command == "pipeline":
        cmd_pipeline(args)
    elif args.command == "suggest-catalog":
        cmd_suggest_catalog(args)
    elif args.command == "validate-sites":
        cmd_validate(args)
    else:
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()
