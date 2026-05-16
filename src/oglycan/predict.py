"""Predictor scaffolding for starter site catalogs.

This module does not implement a state-of-the-art glycosite predictor.
All emitted sites are review-only candidates and are never scored by
core.evaluate().
"""

from __future__ import annotations

import re
import socket
import urllib.error
import urllib.parse
import urllib.request

GLYCOSITE_RESIDUES = ("S", "T")
_PREDICTED_SITE_FIELDS = ("pos", "aa", "p_glycosite", "source")
NETOGLYC_REQUEST_FIELDS = {
    "SEQPASTE": "<sequence placeholder>",
    "outputformat": "short",
}
NETOGLYC_SOURCE = "netoglyc_4.0"


def _normalize_sequence(sequence: str) -> str:
    return re.sub(r"\s+", "", sequence).upper()


def parse_fasta(fasta_text: str) -> list[tuple[str, str]]:
    """Minimal FASTA parser. Returns list of (header, sequence) tuples."""
    records: list[tuple[str, str]] = []
    header: str | None = None
    sequence_lines: list[str] = []

    for raw_line in fasta_text.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append((header, _normalize_sequence("".join(sequence_lines))))
            header = line[1:].strip()
            sequence_lines = []
            continue
        if header is not None:
            sequence_lines.append(line)

    if header is not None:
        records.append((header, _normalize_sequence("".join(sequence_lines))))

    return records


def scan_st_positions(sequence: str) -> list[dict]:
    """Return every Ser/Thr position as a review-only candidate."""
    predicted_sites = []
    for index, aa in enumerate(_normalize_sequence(sequence), start=1):
        if aa in GLYCOSITE_RESIDUES:
            predicted_sites.append(
                {
                    "pos": index,
                    "aa": aa,
                    "p_glycosite": None,
                    "source": "scan_st",
                }
            )
    return predicted_sites


def validate_predicted_sites(sites_list: list[dict]) -> list[str]:
    """Validate optional predicted_sites entries."""
    errors: list[str] = []
    required_fields = set(_PREDICTED_SITE_FIELDS)

    for i, site in enumerate(sites_list):
        site_fields = set(site.keys())
        missing = sorted(required_fields - site_fields)
        extra = sorted(site_fields - required_fields)
        if missing:
            errors.append(
                f"predicted_sites[{i}]: missing field(s): {', '.join(repr(field) for field in missing)}"
            )
        if extra:
            errors.append(
                f"predicted_sites[{i}]: unexpected field(s): {', '.join(repr(field) for field in extra)}"
            )
        if missing or extra:
            continue

        pos = site["pos"]
        source = site["source"]
        if isinstance(pos, bool) or not isinstance(pos, int) or pos < 1:
            errors.append(f"predicted_sites[{i}]: pos must be a positive integer")
        if not isinstance(source, str) or not source:
            errors.append(f"predicted_sites[{i}]: source must be a non-empty string")

        aa = site["aa"]
        if aa not in GLYCOSITE_RESIDUES:
            errors.append(f"predicted_sites[{i}]: aa must be one of {GLYCOSITE_RESIDUES}")

        p_glycosite = site["p_glycosite"]
        if p_glycosite is not None:
            if isinstance(p_glycosite, bool) or not isinstance(p_glycosite, (int, float)):
                errors.append(f"predicted_sites[{i}]: p_glycosite must be a float in [0, 1] or null")
            elif not (0.0 <= float(p_glycosite) <= 1.0):
                errors.append(f"predicted_sites[{i}]: p_glycosite {p_glycosite} out of [0, 1]")

    return errors


def _truncate_netoglyc_text(text: str, *, limit: int = 120) -> str:
    single_line = re.sub(r" +", " ", text.replace("\r", " ").replace("\n", " ").strip())
    if len(single_line) <= limit:
        return single_line
    return f"{single_line[: limit - 3]}..."


def _parse_netoglyc_response(text: str) -> list[dict]:
    """Parse NetOGlyc 4.0 tabular output into predicted_sites rows."""
    predicted_sites: list[dict] = []

    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue

        fields = raw_line.rstrip("\r\n").split("\t")
        if len(fields) < 4:
            raise ValueError(
                "expected at least 4 tab-separated fields, "
                f"got {len(fields)} in line: {_truncate_netoglyc_text(line)}"
            )

        try:
            position = int(fields[0])
        except ValueError as exc:
            raise ValueError(
                f"invalid residue position in line: {_truncate_netoglyc_text(line)}"
            ) from exc

        aa_field = fields[1].strip()
        if not aa_field:
            raise ValueError(
                f"missing residue code in line: {_truncate_netoglyc_text(line)}"
            )
        aa = aa_field[0].upper()

        try:
            score = float(fields[2])
        except ValueError as exc:
            raise ValueError(
                f"invalid glycosite score in line: {_truncate_netoglyc_text(line)}"
            ) from exc

        if position < 1:
            raise ValueError(
                f"invalid residue position in line: {_truncate_netoglyc_text(line)}"
            )
        if not 0.0 <= score <= 1.0:
            raise ValueError(
                f"glycosite score out of [0, 1] in line: {_truncate_netoglyc_text(line)}"
            )
        if aa not in GLYCOSITE_RESIDUES:
            continue

        predicted_sites.append(
            {
                "pos": position,
                "aa": aa,
                "p_glycosite": score,
                "source": NETOGLYC_SOURCE,
            }
        )

    return predicted_sites


def predict_netoglyc(
    sequence: str,
    *,
    threshold: float = 0.5,
    base_url: str = "https://services.healthtech.dtu.dk/services/NetOGlyc-4.0/",
    timeout: float = 120.0,
    opener: urllib.request.OpenerDirector | None = None,
) -> list[dict]:
    """Call DTU's NetOGlyc 4.0 prediction service.

    Returns predicted_sites entries with:
      - pos: 1-indexed position
      - aa: "S" or "T"
      - p_glycosite: float in [0, 1]
      - source: "netoglyc_4.0"

    Only returns sites with p_glycosite >= threshold and aa in {"S", "T"}.
    DTU has historically accepted form-encoded POST fields like
    ``NETOGLYC_REQUEST_FIELDS``; if they revise the service, update those field
    names before first real use and verify against your NetOGlyc deployment.
    """
    normalized_sequence = _normalize_sequence(sequence)
    request_fields = {
        key: normalized_sequence if value == "<sequence placeholder>" else value
        for key, value in NETOGLYC_REQUEST_FIELDS.items()
    }
    payload = urllib.parse.urlencode(request_fields).encode("utf-8")
    request = urllib.request.Request(
        base_url,
        data=payload,
        headers={
            "Accept": "text/plain",
            "Content-Type": "application/x-www-form-urlencoded",
        },
        method="POST",
    )
    open_request = opener.open if opener is not None else urllib.request.urlopen

    try:
        with open_request(request, timeout=timeout) as response:
            status = getattr(response, "status", None)
            if status is None and hasattr(response, "getcode"):
                status = response.getcode()
            body = response.read()
    except urllib.error.HTTPError as exc:
        body_excerpt = _truncate_netoglyc_text(exc.read().decode("utf-8", errors="replace"))
        raise RuntimeError(
            f"NetOGlyc returned HTTP {exc.code}: {body_excerpt or exc.reason}"
        ) from exc
    except urllib.error.URLError as exc:
        reason = exc.reason
        if isinstance(reason, socket.timeout):
            reason = str(reason) or "timed out"
        raise RuntimeError(f"NetOGlyc service unreachable: {reason}") from exc
    except socket.timeout as exc:
        raise RuntimeError(f"NetOGlyc service unreachable: {exc}") from exc

    if status is not None and not 200 <= int(status) < 300:
        body_excerpt = _truncate_netoglyc_text(body.decode("utf-8", errors="replace"))
        raise RuntimeError(f"NetOGlyc returned HTTP {status}: {body_excerpt or '(no body)'}")

    if not body or not body.strip():
        raise RuntimeError("NetOGlyc returned empty response")

    predicted_sites = _parse_netoglyc_response(body.decode("utf-8", errors="replace"))
    if not predicted_sites:
        raise RuntimeError("NetOGlyc response contained no predictions")

    return [
        site for site in predicted_sites if float(site["p_glycosite"]) >= threshold
    ]


def predict_stackoglypred_plm(sequence: str, *args, **kwargs) -> list[dict]:
    """Stub for a future Stack-OglyPred-PLM adapter."""
    raise NotImplementedError(
        "Stack-OglyPred-PLM adapter pending. Requires torch + transformers "
        "(not installed as a runtime dependency). To wire this, load the "
        "published model weights and map per-residue probabilities to "
        "predicted_sites entries with source='stack_oglypred_plm'."
    )


def suggest_catalog(
    sequence: str,
    protein_name: str,
    predictor: str = "scan_st",
    threshold: float = 0.0,
) -> dict:
    """Produce a starter catalog with review-only predicted_sites."""
    normalized_sequence = _normalize_sequence(sequence)
    predictors = {
        "scan_st": scan_st_positions,
        "netoglyc": predict_netoglyc,
        "stackoglypred": predict_stackoglypred_plm,
    }
    if predictor not in predictors:
        raise ValueError(f"unsupported predictor: {predictor}")

    predictor_fn = predictors[predictor]
    if predictor == "netoglyc":
        predictor_kwargs = {"threshold": threshold} if threshold > 0.0 else {}
        predicted_sites = predictor_fn(normalized_sequence, **predictor_kwargs)
    else:
        predicted_sites = predictor_fn(normalized_sequence)
        predicted_sites = [
            site
            for site in predicted_sites
            if site["p_glycosite"] is None or site["p_glycosite"] >= threshold
        ]

    return {
        "glycoprotein": {
            "name": protein_name,
            "description": (
                f"Starter catalog from predict.suggest_catalog ({predictor}). "
                "User must review predicted_sites and promote real observations "
                "into the empty `sites` array."
            ),
            "uniprot": None,
            "reference": None,
        },
        "sites": [],
        "predicted_sites": predicted_sites,
        "n_glycosites": [],
        "localization_threshold": 0.75,
    }
