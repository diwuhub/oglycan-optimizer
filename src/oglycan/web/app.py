"""Scientist-first Streamlit UI for oglycan-optimizer."""

from __future__ import annotations

if __package__ in (None, ""):
    import sys
    from pathlib import Path

    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    __package__ = "oglycan.web"

import copy
import hashlib
import html
import json
import re
import tempfile
from pathlib import Path

import streamlit as st

from ..cli import _catalog_validation_errors, _compare_row, _find_default_params
from ..core import evaluate, load_site_catalog
from ..glycan import glycoform_label, glycoform_mass, render_glycoform_svg
from ..method_card import build_method_card
from ..pilot import load_pilot_canonical, load_pilot_msfragger, validate_pilot
from ..predict import parse_fasta
from ..recalibrate import recalibrate_difficulty, write_catalog
from ..recommend import recommend_next_method
from ..report import json_report
from ..sensitivity import get_nested, set_nested, sensitivity
from ..uncertainty import composite_with_ci

REPO_ROOT = Path(__file__).resolve().parents[3]
SITES_DIR = REPO_ROOT / "sites"
DEFAULT_PARAMS = Path(_find_default_params())
RESOLUTION_OPTIONS = [15000, 30000, 60000, 120000]
PARAM_WIDGET_KEYS = {
    "fragmentation.ethcd_sa_percent": "param_ethcd_sa_percent",
    "fragmentation.collision_energy_nce": "param_collision_energy_nce",
    "ms_acquisition.resolution_ms2": "param_resolution_ms2",
    "lc_gradient.gradient_time_min": "param_gradient_time_min",
    "enzyme_preprocessing.use_sialexo": "param_use_sialexo",
    "enzyme_preprocessing.use_pngasef": "param_use_pngasef",
}
PARAM_HELP = {
    "fragmentation.ethcd_sa_percent": "Supplemental HCD activation during EThcD. Sweet spot 15-25%. Above 25% destroys the c/z radical ions needed for O-glycosite localization. Riley & Coon 2018; JACS Au 2025.",
    "fragmentation.collision_energy_nce": "Normalized collision energy for HCD fragmentation. 28-32 optimal on Orbitrap for oxonium-ion generation. Below 25 gives few diagnostic fragments; above 35 over-fragments larger glycan ions. Reiding et al. 2018.",
    "ms_acquisition.resolution_ms2": "MS2 mass resolving power. Higher resolution = better mass accuracy but longer transient time and slower scan rate. 30K is the sweet spot for glycopeptide LC peaks (10-15 s wide); 60K doubles transient to 15 ms and halves scan rate. Kelstrup et al. 2012.",
    "lc_gradient.gradient_time_min": "LC gradient duration. Peak capacity scales with √(time). 90-120 min typical for EThcD glycopeptide workflows; saturates above ~150 min; below 45 min co-elutes glycoforms. Reiding et al. 2018.",
    "enzyme_preprocessing.use_sialexo": "Sialidase pretreatment. Removes sialic-acid caps that sterically block OpeRATOR access. Yields ~+17% glycopeptide recovery on sialylated substrates. Genovis AN-0042.",
    "enzyme_preprocessing.use_pngasef": "N-glycan-releasing enzyme. Simplifies mixtures that carry both N- and O-glycans (e.g., Etanercept). Reduces spectral complexity and chimera rate.",
}
WHY_PARAMETERS_TEXT = (
    "EThcD method quality depends on balancing radical-driven localization against "
    "oxonium-ion yield: too much supplemental activation or NCE suppresses the c/z "
    "series, while too little leaves glycan diagnostics weak. MS2 resolution improves "
    "mass accuracy, but every jump in resolving power stretches the transient and "
    "costs duty cycle on 10-15 second LC peaks. SialEXO matters upstream because "
    "OpeRATOR access is materially better after terminal sialic acids are removed, so "
    "enzyme pretreatment and fragmentation settings should be tuned together rather "
    "than independently."
)


def _load_json(path: Path | None = None, text: str | None = None) -> dict:
    return json.loads(text) if text is not None else json.loads(path.read_text(encoding="utf-8"))


def _site_table(rows: list[dict]) -> str:
    parts = [
        "<table><thead><tr><th>position</th><th>aa</th><th>core_types</th><th>difficulty</th><th>confidence</th><th>pass</th></tr></thead><tbody>"
    ]
    for row in rows:
        color = "#15803d" if row["pass"] else "#b91c1c"
        parts.append(
            "<tr>"
            f"<td>{row['position']}</td><td>{row['amino_acid']}</td>"
            f"<td>{', '.join(row['core_types'])}</td><td>{row['difficulty']:.2f}</td>"
            f"<td>{row['confidence']:.4f}</td><td style='color:{color};font-weight:600'>{'PASS' if row['pass'] else 'FAIL'}</td>"
            "</tr>"
        )
    parts.append("</tbody></table>")
    return "".join(parts)


def _histogram(values: list[float], bins: int = 16) -> list[dict]:
    low, high = min(values), max(values)
    if low == high:
        return [{"bin": f"{low:.4f}", "count": len(values)}]
    width = (high - low) / bins
    counts = [0] * bins
    for value in values:
        counts[min(int((value - low) / width), bins - 1)] += 1
    return [{"bin": f"{low + i * width:.4f}", "count": counts[i]} for i in range(bins)]


def _editor_rows(data) -> list[dict]:
    if hasattr(data, "to_dict"):
        try:
            return data.to_dict("records")
        except TypeError:
            pass
    if isinstance(data, list):
        return [dict(row) for row in data]
    if isinstance(data, tuple):
        return [dict(row) for row in data]
    return [dict(row) for row in list(data)]


def _is_blank(value) -> bool:
    return value is None or value == "" or value != value


def _split_csv(value) -> list[str]:
    if _is_blank(value):
        return []
    if isinstance(value, list):
        return [str(item).strip() for item in value if str(item).strip()]
    return [item.strip() for item in str(value).split(",") if item.strip()]


def _coerce_editor_pos(value):
    if _is_blank(value):
        return None
    if isinstance(value, bool):
        return value
    if isinstance(value, int):
        return value
    if isinstance(value, float):
        return int(value) if value.is_integer() else value
    text = str(value).strip()
    return int(text) if text.isdigit() else value


def _coerce_editor_difficulty(value):
    if _is_blank(value):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return value


def _site_editor_rows(sites: list[dict], *, use_list_column: bool) -> list[dict]:
    rows = []
    for site in sites or []:
        rows.append(
            {
                "pos": site.get("pos"),
                "aa": site.get("aa", ""),
                "core_types": (
                    list(site.get("core_types", []))
                    if use_list_column
                    else ", ".join(site.get("core_types", []))
                ),
                "difficulty": site.get("difficulty"),
                "domain": site.get("domain", ""),
                "immunogenicity_flags": ", ".join(site.get("immunogenicity_flags", [])),
            }
        )
    return rows


def _normalize_site_editor_rows(rows) -> list[dict]:
    normalized = []
    for row in _editor_rows(rows):
        if all(
            _is_blank(row.get(field))
            for field in (
                "pos",
                "aa",
                "core_types",
                "difficulty",
                "domain",
                "immunogenicity_flags",
            )
        ):
            continue
        normalized.append(
            {
                "pos": _coerce_editor_pos(row.get("pos")),
                "aa": str(row.get("aa", "")).strip().upper(),
                "core_types": _split_csv(row.get("core_types")),
                "difficulty": _coerce_editor_difficulty(row.get("difficulty")),
                "domain": "" if _is_blank(row.get("domain")) else str(row.get("domain")).strip(),
                "immunogenicity_flags": _split_csv(row.get("immunogenicity_flags")),
            }
        )
    return normalized


def _predicted_site_rows(predicted_sites: list[dict]) -> list[dict]:
    rows = []
    for site in predicted_sites or []:
        rows.append(
            {
                "promote": False,
                "pos": site.get("pos"),
                "aa": site.get("aa", ""),
                "p_glycosite": site.get("p_glycosite"),
                "source": site.get("source", ""),
            }
        )
    return rows


def _set_editor_catalog_state(catalog: dict, loaded_name: str) -> None:
    st.session_state["editor_catalog"] = copy.deepcopy(catalog)
    st.session_state["editor_original_catalog"] = copy.deepcopy(catalog)
    st.session_state["editor_loaded_name"] = loaded_name
    st.session_state["editor_source_nonce"] = st.session_state.get("editor_source_nonce", 0) + 1
    st.session_state.pop("editor_validation_errors", None)
    st.session_state.pop("editor_rescore", None)
    st.session_state.pop("editor_notice", None)


def _render_catalog_editor(
    *,
    catalog_name: str,
    catalog: dict,
    params: dict,
) -> None:
    st.session_state.setdefault("editor_catalog", None)
    st.session_state.setdefault("editor_original_catalog", None)
    st.session_state.setdefault("editor_loaded_name", catalog_name)
    st.session_state.setdefault("editor_source_nonce", 0)
    st.session_state.setdefault("editor_uploader_nonce", 0)
    st.session_state.setdefault("editor_catalog_name", None)

    if st.session_state["editor_catalog_name"] != catalog_name:
        st.session_state["editor_catalog_name"] = catalog_name
        st.session_state["editor_uploader_nonce"] += 1
        st.session_state.pop("editor_upload_key", None)
        _set_editor_catalog_state(catalog, catalog_name)

    st.caption("Uses the current catalog unless you upload an override for this editor session.")
    uploaded_catalog = st.file_uploader(
        "Upload catalog JSON",
        type=["json"],
        key=f"editor_catalog_upload_{st.session_state['editor_uploader_nonce']}",
    )
    if uploaded_catalog is not None:
        uploaded_bytes = uploaded_catalog.getvalue()
        upload_key = f"{uploaded_catalog.name}:{hashlib.sha256(uploaded_bytes).hexdigest()}"
        if st.session_state.get("editor_upload_key") != upload_key:
            try:
                uploaded_value = _load_json(text=uploaded_bytes.decode("utf-8"))
            except json.JSONDecodeError as exc:
                st.error(f"Catalog upload failed: {exc}")
            else:
                _set_editor_catalog_state(uploaded_value, uploaded_catalog.name)
                st.session_state["editor_upload_key"] = upload_key

    if st.session_state["editor_catalog"] is None:
        st.session_state["editor_catalog"] = copy.deepcopy(catalog)

    editor_catalog = st.session_state["editor_catalog"]
    st.session_state["editor_catalog"] = editor_catalog
    loaded_name = st.session_state["editor_loaded_name"]

    if notice := st.session_state.pop("editor_notice", None):
        st.success(notice)

    st.markdown(
        "Loaded: "
        f"`{loaded_name}`, "
        f"{len(editor_catalog.get('sites', []))} scored sites, "
        f"{len(editor_catalog.get('predicted_sites') or [])} predicted candidates"
    )

    st.subheader("Edit Scored `sites[]`")
    use_list_column = hasattr(st.column_config, "ListColumn")
    core_types_column = (
        st.column_config.ListColumn("core_types", help="Allowed core glycans")
        if use_list_column
        else st.column_config.TextColumn("core_types", help="Comma-separated core glycans")
    )
    edited_sites = st.data_editor(
        _site_editor_rows(editor_catalog.get("sites", []), use_list_column=use_list_column),
        num_rows="dynamic",
        use_container_width=True,
        key=f"editor_sites_table_{st.session_state['editor_source_nonce']}",
        column_config={
            "pos": st.column_config.NumberColumn("pos", min_value=1, step=1, format="%d"),
            "aa": st.column_config.SelectboxColumn("aa", options=["S", "T"]),
            "core_types": core_types_column,
            "difficulty": st.column_config.NumberColumn(
                "difficulty",
                min_value=0.0,
                max_value=1.0,
                step=0.05,
                format="%.2f",
            ),
            "domain": st.column_config.TextColumn("domain"),
            "immunogenicity_flags": st.column_config.TextColumn(
                "immunogenicity_flags",
                help="Comma-separated flags",
            ),
        },
    )
    normalized_sites = _normalize_site_editor_rows(edited_sites)
    if normalized_sites != editor_catalog.get("sites", []):
        editor_catalog["sites"] = normalized_sites
        st.session_state.pop("editor_validation_errors", None)
        st.session_state.pop("editor_rescore", None)

    st.subheader("Promote `predicted_sites` -> `sites`")
    predicted_sites = editor_catalog.get("predicted_sites") or []
    if predicted_sites:
        predicted_rows = st.data_editor(
            _predicted_site_rows(predicted_sites),
            hide_index=True,
            use_container_width=True,
            disabled=["pos", "aa", "p_glycosite", "source"],
            key=f"editor_predicted_table_{st.session_state['editor_source_nonce']}",
            column_config={
                "promote": st.column_config.CheckboxColumn("Promote"),
                "pos": st.column_config.NumberColumn("pos", format="%d"),
                "aa": st.column_config.TextColumn("aa"),
                "p_glycosite": st.column_config.NumberColumn("p_glycosite", format="%.3f"),
                "source": st.column_config.TextColumn("source"),
            },
        )
        if st.button("Promote selected"):
            selected_rows = [row for row in _editor_rows(predicted_rows) if row.get("promote")]
            if not selected_rows:
                st.warning("Select at least one predicted site to promote.")
            else:
                selected_positions = {
                    (_coerce_editor_pos(row.get("pos")), str(row.get("aa", "")).strip().upper())
                    for row in selected_rows
                }
                editor_catalog["sites"] = editor_catalog.get("sites", []) + [
                    {
                        "pos": _coerce_editor_pos(row.get("pos")),
                        "aa": str(row.get("aa", "")).strip().upper(),
                        "core_types": ["Core1"],
                        "difficulty": 0.5,
                        "domain": "",
                        "immunogenicity_flags": [],
                    }
                    for row in selected_rows
                ]
                editor_catalog["predicted_sites"] = [
                    site
                    for site in predicted_sites
                    if (site.get("pos"), str(site.get("aa", "")).strip().upper())
                    not in selected_positions
                ]
                st.session_state["editor_notice"] = f"Promoted {len(selected_rows)} predicted site(s)."
                st.session_state["editor_source_nonce"] += 1
                st.session_state.pop("editor_validation_errors", None)
                st.session_state.pop("editor_rescore", None)
                st.rerun()
    else:
        st.info("No predicted candidates remain.")

    st.subheader("Validate")
    if st.button("Validate"):
        st.session_state["editor_validation_errors"] = _catalog_validation_errors(editor_catalog)
    if "editor_validation_errors" in st.session_state:
        errors = st.session_state["editor_validation_errors"]
        if errors:
            for error in errors:
                st.error(error)
        else:
            st.success("All sites valid.")

    st.subheader("Save")
    download_name = f"{Path(loaded_name).stem}_edited.json"
    st.download_button(
        "Download catalog",
        data=json.dumps(editor_catalog, indent=2) + "\n",
        file_name=download_name,
        mime="application/json",
    )
    if st.button("Re-score with current params"):
        errors = _catalog_validation_errors(editor_catalog)
        st.session_state["editor_validation_errors"] = errors
        if not errors:
            st.session_state["editor_rescore"] = {
                "original": evaluate(st.session_state["editor_original_catalog"], params),
                "edited": evaluate(editor_catalog, params),
            }
        else:
            st.session_state.pop("editor_rescore", None)
    if "editor_rescore" in st.session_state:
        original_col, edited_col, delta_col = st.columns(3)
        original_score = st.session_state["editor_rescore"]["original"]["composite_score"]
        edited_score = st.session_state["editor_rescore"]["edited"]["composite_score"]
        original_col.metric("Original composite", f"{original_score:.4f}")
        edited_col.metric("Edited composite", f"{edited_score:.4f}")
        delta_col.metric("Delta", f"{edited_score - original_score:.4f}")


def _normalize_name(text: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", text.lower())


def _normalize_sequence_text(text: str) -> str:
    stripped = (text or "").strip()
    if not stripped:
        return ""
    if stripped.startswith(">"):
        records = parse_fasta(stripped)
        if records:
            return records[0][1]
    return re.sub(r"\s+", "", stripped).upper()


def _display_name_from_header(header: str) -> str:
    header = (header or "").strip()
    if not header:
        return ""
    if "|" in header:
        fields = [field.strip() for field in header.split("|") if field.strip()]
        if fields:
            header = fields[-1]
    return header.split()[0] if header.split() else header


def _catalog_entries(catalog_paths: dict[str, Path]) -> list[dict]:
    entries = []
    for catalog_name, path in sorted(catalog_paths.items()):
        catalog = load_site_catalog(str(path))
        protein_name = catalog.get("glycoprotein", {}).get("name", Path(catalog_name).stem)
        aliases = {
            _normalize_name(protein_name),
            _normalize_name(Path(catalog_name).stem),
        }
        entries.append(
            {
                "catalog_name": catalog_name,
                "catalog": catalog,
                "path": path,
                "protein_name": protein_name,
                "aliases": {alias for alias in aliases if alias},
            }
        )
    return entries


def _find_matching_catalog_name(protein_name: str, entries: list[dict]) -> str | None:
    query = _normalize_name(protein_name)
    if not query:
        return None

    best_name = None
    best_score = None
    for entry in entries:
        for alias in entry["aliases"]:
            if query == alias:
                score = (0, -len(alias))
            elif alias in query:
                score = (1, -len(alias))
            elif query in alias:
                score = (2, -len(alias))
            else:
                continue
            if best_score is None or score < best_score:
                best_score = score
                best_name = entry["catalog_name"]
    return best_name


def _site_result_lookup(site_results: list[dict]) -> dict[tuple[int, str], dict]:
    return {(row["position"], row["amino_acid"].upper()): row for row in site_results}


def _site_score_value(pos: int, aa: str, site_scores: dict) -> float | None:
    for key in ((pos, aa.upper()), pos):
        if key not in site_scores:
            continue
        value = site_scores[key]
        if isinstance(value, dict):
            value = value.get("confidence")
        if value is None:
            return None
        return float(value)
    return None


def _site_lookup(catalog: dict) -> set[tuple[int, str]]:
    return {
        (int(site["pos"]), str(site["aa"]).strip().upper())
        for site in catalog.get("sites", [])
    }


def _predicted_lookup(catalog: dict) -> set[tuple[int, str]]:
    return {
        (int(site["pos"]), str(site["aa"]).strip().upper())
        for site in catalog.get("predicted_sites") or []
    }


def _classify_residue(pos: int, aa: str, catalog: dict, site_scores: dict) -> str:
    """Return one of: 'none' (non-S/T), 'unmarked' (S/T no evidence),
    'predicted' (in predicted_sites), 'pass' (in sites, confidence >= threshold),
    'fail' (in sites, confidence < threshold)."""
    aa = aa.upper()
    if aa not in {"S", "T"}:
        return "none"

    key = (pos, aa)
    threshold = float(catalog.get("localization_threshold", 0.75))
    if key in _site_lookup(catalog):
        confidence = _site_score_value(pos, aa, site_scores)
        return "pass" if confidence is not None and confidence >= threshold else "fail"
    if key in _predicted_lookup(catalog):
        return "predicted"
    return "unmarked"


def render_annotated_sequence_html(sequence: str, catalog: dict, site_scores: dict) -> str:
    """Return an HTML string ready for st.markdown(unsafe_allow_html=True)."""
    normalized = _normalize_sequence_text(sequence)
    if not normalized:
        return (
            "<div style='font-family: monospace; color: #6b7280;'>"
            "Upload a FASTA or paste a sequence to render the annotated viewer."
            "</div>"
        )

    style = """
<style>
.oglycan-sequence { font-family: monospace; color: #111827; }
.oglycan-line { margin-bottom: 0.45rem; }
.oglycan-ruler { color: #9ca3af; line-height: 1.1; }
.oglycan-ruler-block { display: inline-block; text-align: right; }
.oglycan-residue { display: inline-block; width: 1ch; text-align: center; }
.oglycan-none { color: #111827; }
.oglycan-unmarked { background: #e5e7eb; color: #4b5563; font-style: italic; border-radius: 2px; }
.oglycan-predicted { background: #fde68a; color: #111827; border-radius: 2px; }
.oglycan-pass { background: #86efac; color: #14532d; border-radius: 2px; }
.oglycan-fail { background: #fca5a5; color: #7f1d1d; border-radius: 2px; }
.oglycan-legend { margin-top: 0.65rem; color: #374151; }
.oglycan-legend-item { display: inline-flex; align-items: center; margin-right: 1rem; margin-bottom: 0.25rem; }
.oglycan-legend-swatch { width: 0.9em; height: 0.9em; border: 1px solid #4b5563; margin-right: 0.35rem; display: inline-block; }
</style>
"""

    parts = [style, "<div class='oglycan-sequence'>"]
    for line_start in range(0, len(normalized), 80):
        line = normalized[line_start:line_start + 80]
        parts.append("<div class='oglycan-line'>")
        parts.append("<div class='oglycan-ruler'>")
        for block_start in range(0, len(line), 10):
            block_len = min(10, len(line) - block_start)
            line_end = line_start + block_start + block_len
            parts.append(
                f"<span class='oglycan-ruler-block' style='width: {block_len}ch'>{line_end}</span>"
            )
        parts.append("</div><div>")
        for offset, aa in enumerate(line, start=line_start + 1):
            status = _classify_residue(offset, aa, catalog, site_scores)
            parts.append(
                f"<span class='oglycan-residue oglycan-{status}'>{html.escape(aa)}</span>"
            )
        parts.append("</div></div>")
    parts.extend(
        [
            "<div class='oglycan-legend'>",
            "<span class='oglycan-legend-item'><span class='oglycan-legend-swatch oglycan-unmarked'></span>Ser/Thr with no evidence</span>",
            "<span class='oglycan-legend-item'><span class='oglycan-legend-swatch oglycan-predicted'></span>Predicted candidate</span>",
            "<span class='oglycan-legend-item'><span class='oglycan-legend-swatch oglycan-pass'></span>Confirmed, passes threshold</span>",
            "<span class='oglycan-legend-item'><span class='oglycan-legend-swatch oglycan-fail'></span>Confirmed, below threshold</span>",
            "</div></div>",
        ]
    )
    return "".join(parts)


def _site_option_label(site: dict) -> str:
    domain = str(site.get("domain", "")).strip()
    base = f"{site['aa']}{site['pos']}"
    return f"{base} ({domain})" if domain else base


def _site_detail_table_html(site: dict, site_result: dict | None, threshold: float) -> str:
    confidence = None if site_result is None else float(site_result["confidence"])
    passed = False if site_result is None else bool(site_result["pass"])
    parts = [
        "<table><thead><tr><th>SNFG symbol</th><th>Label</th><th>Mass (Da)</th><th>Confidence</th></tr></thead><tbody>"
    ]
    for core_type in site.get("core_types", []):
        symbol = render_glycoform_svg(core_type, symbol_px=18, spacing_px=6)
        label = glycoform_label(core_type)
        mass = glycoform_mass(core_type)
        if confidence is None:
            confidence_html = "n/a"
        else:
            mark = " &#10003;" if passed and confidence >= threshold else ""
            confidence_html = f"{confidence:.4f}{mark}"
        parts.append(
            "<tr>"
            f"<td>{symbol}</td>"
            f"<td>{html.escape(label)}</td>"
            f"<td>{mass:.4f}</td>"
            f"<td>{confidence_html}</td>"
            "</tr>"
        )
    parts.append("</tbody></table>")
    return "".join(parts)


def _weakest_site(evaluation: dict) -> dict | None:
    rows = evaluation.get("site_results", [])
    return min(rows, key=lambda item: item["confidence"]) if rows else None


def _confidence_trend_text(before_eval: dict, after_eval: dict) -> str:
    before_lookup = _site_result_lookup(before_eval.get("site_results", []))
    improved = 0
    worsened = 0
    for key, after_row in _site_result_lookup(after_eval.get("site_results", [])).items():
        before_row = before_lookup.get(key)
        if before_row is None:
            continue
        delta = float(after_row["confidence"]) - float(before_row["confidence"])
        if delta > 0.0001:
            improved += 1
        elif delta < -0.0001:
            worsened += 1
    if improved > worsened:
        return "Confidence improves for most sites."
    if worsened > improved:
        return "Confidence drops slightly for most sites."
    return "Per-site confidence stays broadly similar."


def _recommendation_blurb(
    *,
    before_eval: dict,
    after_eval: dict,
    recommendation: dict,
    constraints: dict,
) -> str:
    changes = recommendation.get("changes") or []
    if not changes:
        return (
            f"The current method already sits at the best composite for this search. "
            f"{after_eval['sites_localized']} of {after_eval['total_sites']} sites pass."
        )

    primary_change = changes[0]
    path = primary_change["path"].split(".")[-1]
    prefix = "With no hard runtime or composite constraint, "
    if constraints.get("max_runtime_min"):
        prefix = f"Under your {constraints['max_runtime_min']:.0f}-minute runtime cap, "
    elif constraints.get("min_composite"):
        prefix = f"To target a minimum engineering composite of {constraints['min_composite']:.2f}, "
    return (
        f"{prefix}the best method switches {path} from {primary_change['before']} to "
        f"{primary_change['after']}. {after_eval['sites_localized']} of "
        f"{after_eval['total_sites']} sites still pass (was {before_eval['sites_localized']}). "
        f"{_confidence_trend_text(before_eval, after_eval)}"
    )


def _pilot_summary_rows(pilot_result: dict) -> list[dict]:
    rows = []
    for site in sorted(pilot_result.get("sites", []), key=lambda item: (item["pos"], item["aa"])):
        rows.append(
            {
                "site": f"{site['aa']}{site['pos']}",
                "localization": site["observed_localization"],
                "n_spectra": site["n_spectra"],
                "glycoforms": ", ".join(site["observed_glycoforms"]),
            }
        )
    return rows


def _sequence_warnings(sequence: str, catalog: dict) -> list[str]:
    warnings = []
    if not sequence:
        return warnings
    max_pos = max((site["pos"] for site in catalog.get("sites", [])), default=0)
    if max_pos > len(sequence):
        warnings.append(
            f"Sequence length {len(sequence)} does not reach the highest catalog site ({max_pos})."
        )
    for site in catalog.get("sites", []):
        pos = int(site["pos"])
        if pos > len(sequence):
            continue
        observed = sequence[pos - 1].upper()
        expected = str(site["aa"]).upper()
        if observed != expected:
            warnings.append(
                f"Catalog site {expected}{pos} disagrees with the input sequence residue ({observed})."
            )
    return warnings


def _ensure_session_state(default_catalog_name: str, default_catalog: dict) -> None:
    if "catalog_name" not in st.session_state:
        st.session_state["catalog_name"] = default_catalog_name
    if "params" not in st.session_state:
        _set_params_state(_load_json(DEFAULT_PARAMS))
    else:
        for path, key in PARAM_WIDGET_KEYS.items():
            if key not in st.session_state:
                st.session_state[key] = get_nested(st.session_state["params"], path)
    st.session_state.setdefault(
        "protein_name",
        default_catalog.get("glycoprotein", {}).get("name", Path(default_catalog_name).stem),
    )
    st.session_state.setdefault("sequence_text", "")
    st.session_state.setdefault("show_recommendation_panel", False)
    st.session_state.setdefault("recommend_max_runtime_min", 0.0)
    st.session_state.setdefault("recommend_min_composite", 0.0)


def _set_params_state(params: dict) -> None:
    st.session_state["params"] = copy.deepcopy(params)
    for path, key in PARAM_WIDGET_KEYS.items():
        st.session_state[key] = get_nested(params, path)


def _apply_pending_state() -> None:
    pending_catalog_name = st.session_state.pop("pending_catalog_name", None)
    if pending_catalog_name is not None:
        st.session_state["catalog_name"] = pending_catalog_name

    pending_params = st.session_state.pop("pending_params", None)
    if pending_params is not None:
        _set_params_state(pending_params)


def _current_params() -> dict:
    params = copy.deepcopy(st.session_state["params"])
    for path, key in PARAM_WIDGET_KEYS.items():
        set_nested(params, path, st.session_state[key])
    st.session_state["params"] = copy.deepcopy(params)
    return params


def main() -> None:
    st.set_page_config(page_title="oglycan-optimizer", layout="wide")
    catalog_paths = {path.name: path for path in sorted(SITES_DIR.glob("*.json"))}
    if not catalog_paths:
        st.error("No catalogs found in sites/.")
        return

    default_catalog_name = "etanercept.json" if "etanercept.json" in catalog_paths else next(iter(catalog_paths))
    default_catalog = load_site_catalog(str(catalog_paths[default_catalog_name]))
    _ensure_session_state(default_catalog_name, default_catalog)
    _apply_pending_state()

    current_catalog_name = st.session_state["catalog_name"]
    if current_catalog_name not in catalog_paths:
        current_catalog_name = default_catalog_name
        st.session_state["catalog_name"] = default_catalog_name
    catalog_path = catalog_paths[current_catalog_name]
    catalog = load_site_catalog(str(catalog_path))
    catalog_entries = _catalog_entries(catalog_paths)
    params = _current_params()
    evaluation = evaluate(catalog, params)
    site_scores = _site_result_lookup(evaluation["site_results"])
    method_card = build_method_card(catalog, params, evaluation)
    report_json_text = json_report(
        evaluation,
        source_catalog=str(catalog_path),
        full_params=params,
    )

    source_key = f"{current_catalog_name}:{json.dumps(params, sort_keys=True)}"
    if st.session_state.get("recommendation_source_key") != source_key:
        st.session_state.pop("recommendation_payload", None)

    analyze_tab, advanced_tab = st.tabs(["Analyze", "Advanced"])

    with analyze_tab:
        st.title("oglycan-optimizer")
        st.caption(
            f"{evaluation['molecule']} with {evaluation['total_sites']} cataloged O-glycosites"
        )

        st.header("1. Protein input")
        uploaded_fasta = st.file_uploader("FASTA file uploader", type=["fa", "fasta", "txt"])
        if uploaded_fasta is not None:
            uploaded_bytes = uploaded_fasta.getvalue()
            upload_sig = f"{uploaded_fasta.name}:{hashlib.sha256(uploaded_bytes).hexdigest()}"
            if st.session_state.get("fasta_upload_sig") != upload_sig:
                records = parse_fasta(uploaded_bytes.decode("utf-8", errors="replace"))
                if not records:
                    st.error("No FASTA records found in the uploaded file.")
                else:
                    header, uploaded_sequence = records[0]
                    st.session_state["sequence_text"] = uploaded_sequence
                    st.session_state["protein_name"] = _display_name_from_header(header)
                    st.session_state["fasta_upload_sig"] = upload_sig
                    if len(records) > 1:
                        st.info(
                            f"Loaded the first FASTA record from {len(records)} records in the file."
                        )

        sequence_text = st.text_area(
            "Raw sequence text area",
            key="sequence_text",
            height=140,
            help="Paste a raw protein sequence here if you do not have a FASTA file.",
        )
        input_cols = st.columns([1.4, 1.1, 0.9])
        with input_cols[0]:
            protein_name = st.text_input("Protein name field", key="protein_name")
        with input_cols[1]:
            st.selectbox("Current catalog", list(catalog_paths), key="catalog_name")
        with input_cols[2]:
            st.write("")
            if st.button("Load matching catalog"):
                match = _find_matching_catalog_name(protein_name, catalog_entries)
                if match is None:
                    st.session_state["catalog_match_notice"] = (
                        "warning",
                        f"No catalog match found for {protein_name or '(blank name)'}",
                    )
                else:
                    st.session_state["catalog_match_notice"] = (
                        "success",
                        f"Matched {protein_name} to {match}",
                    )
                    st.session_state["pending_catalog_name"] = match
                st.rerun()

        if notice := st.session_state.pop("catalog_match_notice", None):
            level, message = notice
            if level == "success":
                st.success(message)
            else:
                st.warning(message)
        st.caption(f"Using catalog `{current_catalog_name}` for site scoring.")

        sequence = _normalize_sequence_text(sequence_text)

        st.header("2. Sequence viewer")
        warnings = _sequence_warnings(sequence, catalog)
        for warning in warnings:
            st.warning(warning)
        st.markdown(
            render_annotated_sequence_html(sequence, catalog, site_scores),
            unsafe_allow_html=True,
        )

        st.header("3. Inspect a site")
        confirmed_sites = catalog.get("sites", [])
        if confirmed_sites:
            selected_label = st.selectbox(
                "Confirmed site",
                [_site_option_label(site) for site in confirmed_sites],
                key=f"site_select_{current_catalog_name}",
            )
            selected_site = next(
                site for site in confirmed_sites if _site_option_label(site) == selected_label
            )
            selected_result = site_scores.get((selected_site["pos"], selected_site["aa"].upper()))
            domain = str(selected_site.get("domain", "")).strip() or "unassigned"
            flags = ", ".join(selected_site.get("immunogenicity_flags", [])) or "none"
            st.markdown(f"**Site {selected_site['aa']}{selected_site['pos']} ({html.escape(domain)} domain)**")
            st.markdown(
                f"Difficulty: {float(selected_site['difficulty']):.2f}  "
                f"Domain: {html.escape(domain)}  "
                f"Flags: {html.escape(flags)}"
            )
            st.markdown(
                _site_detail_table_html(
                    selected_site,
                    selected_result,
                    float(catalog.get("localization_threshold", 0.75)),
                ),
                unsafe_allow_html=True,
            )
            st.caption(
                "Per-glycoform confidence currently aggregates to the site level in this version."
            )
        else:
            st.info("This catalog has no confirmed sites yet.")

        st.header("4. Acquisition parameters")
        param_left, param_right = st.columns(2)
        with param_left:
            st.slider(
                "EThcD SA %  (supplemental HCD activation)",
                min_value=0,
                max_value=50,
                step=1,
                value=int(st.session_state[PARAM_WIDGET_KEYS["fragmentation.ethcd_sa_percent"]]),
                key=PARAM_WIDGET_KEYS["fragmentation.ethcd_sa_percent"],
                help=PARAM_HELP["fragmentation.ethcd_sa_percent"],
            )
            st.select_slider(
                "MS2 resolution",
                options=RESOLUTION_OPTIONS,
                value=int(st.session_state[PARAM_WIDGET_KEYS["ms_acquisition.resolution_ms2"]]),
                key=PARAM_WIDGET_KEYS["ms_acquisition.resolution_ms2"],
                help=PARAM_HELP["ms_acquisition.resolution_ms2"],
            )
            st.checkbox(
                "SialEXO pretreatment",
                value=bool(st.session_state[PARAM_WIDGET_KEYS["enzyme_preprocessing.use_sialexo"]]),
                key=PARAM_WIDGET_KEYS["enzyme_preprocessing.use_sialexo"],
                help=PARAM_HELP["enzyme_preprocessing.use_sialexo"],
            )
        with param_right:
            st.slider(
                "HCD NCE  (normalized collision energy)",
                min_value=10,
                max_value=50,
                step=1,
                value=int(st.session_state[PARAM_WIDGET_KEYS["fragmentation.collision_energy_nce"]]),
                key=PARAM_WIDGET_KEYS["fragmentation.collision_energy_nce"],
                help=PARAM_HELP["fragmentation.collision_energy_nce"],
            )
            st.slider(
                "LC gradient (min)",
                min_value=30,
                max_value=240,
                step=15,
                value=int(st.session_state[PARAM_WIDGET_KEYS["lc_gradient.gradient_time_min"]]),
                key=PARAM_WIDGET_KEYS["lc_gradient.gradient_time_min"],
                help=PARAM_HELP["lc_gradient.gradient_time_min"],
            )
            st.checkbox(
                "PNGase F pretreatment",
                value=bool(st.session_state[PARAM_WIDGET_KEYS["enzyme_preprocessing.use_pngasef"]]),
                key=PARAM_WIDGET_KEYS["enzyme_preprocessing.use_pngasef"],
                help=PARAM_HELP["enzyme_preprocessing.use_pngasef"],
            )
        with st.expander("Why these parameters?"):
            st.write(WHY_PARAMETERS_TEXT)

        st.header("5. Result summary")
        if evaluation["total_sites"]:
            st.write(
                f"{evaluation['sites_localized']} of {evaluation['total_sites']} sites pass "
                f"at {evaluation['localization_threshold']:.2f} threshold."
            )
            weakest_site = _weakest_site(evaluation)
            if weakest_site is not None:
                st.write(
                    f"Weakest: {weakest_site['amino_acid']}{weakest_site['position']} "
                    f"at {weakest_site['confidence']:.4f}."
                )
        else:
            st.write("No confirmed sites are available in this catalog yet.")
        st.markdown(
            f"<div style='color:#6b7280; font-size:0.9rem'>(engineering composite: {evaluation['composite_score']:.4f})</div>",
            unsafe_allow_html=True,
        )

        st.header("6. Actions")
        if st.button("Recommend better method"):
            st.session_state["show_recommendation_panel"] = True
        st.caption(
            "Raw JSON / markdown exports live in the **Advanced** tab; this view stays focused on the scientific result."
        )

        if st.session_state.get("show_recommendation_panel"):
            st.subheader("Method recommendation")
            rec_cols = st.columns(2)
            with rec_cols[0]:
                max_runtime = st.number_input(
                    "max_runtime_min",
                    min_value=0.0,
                    step=15.0,
                    key="recommend_max_runtime_min",
                )
            with rec_cols[1]:
                min_composite = st.number_input(
                    "min_composite",
                    min_value=0.0,
                    max_value=1.0,
                    step=0.01,
                    key="recommend_min_composite",
                )

            if st.button("Generate recommendation"):
                constraints = {}
                if max_runtime > 0:
                    constraints["max_runtime_min"] = max_runtime
                if min_composite > 0:
                    constraints["min_composite"] = min_composite
                recommendation = recommend_next_method(catalog, params, constraints or None)
                recommended_eval = evaluate(catalog, recommendation["recommended_params"])
                st.session_state["recommendation_payload"] = {
                    "recommendation": recommendation,
                    "recommended_evaluation": recommended_eval,
                    "constraints": constraints,
                }
                st.session_state["recommendation_source_key"] = source_key

            payload = st.session_state.get("recommendation_payload")
            if payload is not None:
                recommendation = payload["recommendation"]
                recommended_eval = payload["recommended_evaluation"]
                constraints = payload["constraints"]
                st.write(
                    _recommendation_blurb(
                        before_eval=evaluation,
                        after_eval=recommended_eval,
                        recommendation=recommendation,
                        constraints=constraints,
                    )
                )
                changes = recommendation["changes"] or [
                    {"path": "(none)", "before": "", "after": ""}
                ]
                st.dataframe(changes, use_container_width=True)
                for violation in recommendation["constraint_violations"]:
                    st.warning(violation)
                if st.button("Apply recommended params"):
                    st.session_state["pending_params"] = recommendation["recommended_params"]
                    st.rerun()

    with advanced_tab:
        st.subheader("Engineering Score Details")
        score_left, score_right = st.columns([1, 2])
        score_left.metric("Composite score", f"{evaluation['composite_score']:.4f}")
        score_right.metric(
            "Sites localized",
            f"{evaluation['sites_localized']}/{evaluation['total_sites']}",
        )
        st.markdown(_site_table(evaluation["site_results"]), unsafe_allow_html=True)

        st.divider()
        st.subheader("Uncertainty")
        ci = composite_with_ci(catalog, params, n_samples=500)
        left, middle, right = st.columns(3)
        left.metric("Point", f"{ci['point']:.4f}")
        middle.metric("CI low", f"{ci['ci_low']:.4f}")
        right.metric("CI high", f"{ci['ci_high']:.4f}")
        st.vega_lite_chart(
            {
                "data": {"values": _histogram(ci["samples"])},
                "mark": "bar",
                "encoding": {
                    "x": {"field": "bin", "type": "ordinal", "sort": None, "title": "sample bins"},
                    "y": {"field": "count", "type": "quantitative"},
                },
            },
            use_container_width=True,
        )

        st.divider()
        st.subheader("Sensitivity")
        st.dataframe(sensitivity(catalog, params)[:10], use_container_width=True)

        st.divider()
        st.subheader("Cross-Protein Compare")
        compare_catalogs = st.multiselect(
            "Catalogs to compare",
            list(catalog_paths),
            default=list(catalog_paths),
            key="advanced_compare_catalogs",
        )
        include_compare_sensitivity = st.checkbox(
            "Include top sensitivity movers",
            key="advanced_compare_include_sensitivity",
        )
        if compare_catalogs:
            compare_rows = []
            for compare_name in compare_catalogs:
                compare_catalog = load_site_catalog(str(catalog_paths[compare_name]))
                compare_eval = evaluate(compare_catalog, params)
                top_rows = sensitivity(compare_catalog, params)[:3] if include_compare_sensitivity else None
                compare_rows.append(
                    _compare_row(
                        str(catalog_paths[compare_name]),
                        compare_catalog,
                        compare_eval,
                        top_sensitivity_rows=top_rows,
                    )
                )
            st.dataframe(compare_rows, use_container_width=True)
        else:
            st.info("Select at least one catalog to compare.")

        st.divider()
        st.subheader("Catalog Editor")
        _render_catalog_editor(catalog_name=current_catalog_name, catalog=catalog, params=params)

        st.divider()
        st.subheader("Pilot Ingestion")
        uploaded_pilot = st.file_uploader(
            "Upload pilot evidence (JSON or TSV)",
            type=["json", "tsv"],
            key="advanced_pilot_upload",
        )
        pilot_result = None
        if uploaded_pilot is not None:
            suffix = Path(uploaded_pilot.name).suffix.lower()
            with tempfile.NamedTemporaryFile(delete=False, suffix=suffix) as handle:
                handle.write(uploaded_pilot.getvalue())
                temp_path = Path(handle.name)
            try:
                if suffix == ".json":
                    pilot_result = load_pilot_canonical(str(temp_path))
                else:
                    pilot_result = load_pilot_msfragger(
                        str(temp_path),
                        glycoprotein_name=catalog["glycoprotein"]["name"],
                        sequence=sequence or None,
                    )
                errors = validate_pilot(pilot_result)
                if errors:
                    raise ValueError("; ".join(errors))
            except Exception as exc:
                st.error(f"Pilot ingestion failed: {exc}")
            finally:
                temp_path.unlink(missing_ok=True)

        if pilot_result is not None:
            st.success(
                f"Pilot loaded: {pilot_result['source']} ({len(pilot_result['sites'])} sites)"
            )
            st.dataframe(_pilot_summary_rows(pilot_result), use_container_width=True)
            if st.button("Recalibrate this catalog", key=f"recalibrate_{current_catalog_name}"):
                recalibrated = recalibrate_difficulty(catalog, pilot_result)
                with tempfile.NamedTemporaryFile(
                    delete=False,
                    prefix=f"{Path(current_catalog_name).stem}_recalibrated_",
                    suffix=".json",
                ) as handle:
                    output_path = Path(handle.name)
                write_catalog(recalibrated, str(output_path))
                st.session_state["pilot_recalibration"] = {
                    "catalog_name": current_catalog_name,
                    "path": str(output_path),
                    "json": output_path.read_text(encoding="utf-8"),
                }

        recalibration = st.session_state.get("pilot_recalibration")
        if recalibration and recalibration.get("catalog_name") == current_catalog_name:
            st.success(f"Recalibrated catalog written to {recalibration['path']}")
            st.download_button(
                "Download recalibrated catalog",
                data=recalibration["json"],
                file_name=f"{Path(current_catalog_name).stem}_recalibrated.json",
                mime="application/json",
            )

        st.divider()
        st.subheader("JSON Download")
        uploaded_params = st.file_uploader(
            "Upload acquisition params JSON",
            type=["json"],
            key="advanced_params_upload",
        )
        if uploaded_params is not None and st.button("Load uploaded params"):
            try:
                loaded_params = _load_json(text=uploaded_params.getvalue().decode("utf-8"))
            except json.JSONDecodeError as exc:
                st.error(f"Params upload failed: {exc}")
            else:
                st.session_state["pending_params"] = loaded_params
                st.rerun()

        download_cols = st.columns(3)
        with download_cols[0]:
            st.download_button(
                "Download params JSON",
                data=json.dumps(params, indent=2) + "\n",
                file_name="acquisition_params.json",
                mime="application/json",
            )
        with download_cols[1]:
            st.download_button(
                "Download report JSON",
                data=report_json_text + "\n",
                file_name="report.json",
                mime="application/json",
            )
        with download_cols[2]:
            st.download_button(
                "Download method card (Markdown)",
                data=method_card,
                file_name=f"{Path(current_catalog_name).stem}_method_card.md",
                mime="text/markdown",
            )


if __name__ == "__main__":
    main()
