"""Glycan chemistry: composition, monoisotopic mass, and SNFG rendering.

The glycoforms named in site catalogs ("Core1", "Core1_Sia", ...) resolve to
minimal canonical compositions for their structural class. Real glycopeptide
MS data may show additional species (longer extensions, multi-antenna, etc.);
the tool presents these as scaffolds for user review, not as ground truth.
"""

from __future__ import annotations

MONOISOTOPIC_RESIDUE_MASS = {
    "GalNAc": 203.079372,
    "GlcNAc": 203.079372,
    "Gal": 162.052824,
    "Man": 162.052824,
    "NeuAc": 291.095416,
    "Fuc": 146.057909,
}

CORE_TYPE_COMPOSITION = {
    "Core1_GalNAc": {"GalNAc": 1},
    "Core1": {"GalNAc": 1, "Gal": 1},
    "Core1_Sia": {"GalNAc": 1, "Gal": 1, "NeuAc": 1},
    "Core2": {"GalNAc": 1, "Gal": 1, "GlcNAc": 1},
    "Core2_Fuc": {"GalNAc": 1, "Gal": 1, "GlcNAc": 1, "Fuc": 1},
    "Core3": {"GalNAc": 1, "GlcNAc": 1},
}

CORE_TYPE_LABEL = {
    "Core1_GalNAc": "Tn (GalNAc)",
    "Core1": "T (Gal-GalNAc)",
    "Core1_Sia": "sialyl-T",
    "Core2": "Core2",
    "Core2_Fuc": "Core2 + Fuc",
    "Core3": "Core3",
}

SNFG_SPEC = {
    "GalNAc": {"shape": "square", "fill": "#F4C92F"},
    "GlcNAc": {"shape": "square", "fill": "#2E77B5"},
    "Gal": {"shape": "circle", "fill": "#F4C92F"},
    "Man": {"shape": "circle", "fill": "#00A651"},
    "NeuAc": {"shape": "diamond", "fill": "#A6519F"},
    "Fuc": {"shape": "triangle", "fill": "#E91313"},
}

_LAYOUTS = {
    "Core1_GalNAc": {
        "nodes": [{"residue": "GalNAc", "x": 0, "y": 0}],
        "edges": [],
    },
    "Core1": {
        "nodes": [
            {"residue": "Gal", "x": 0, "y": 0},
            {"residue": "GalNAc", "x": 1, "y": 0},
        ],
        "edges": [(0, 1)],
    },
    "Core1_Sia": {
        "nodes": [
            {"residue": "NeuAc", "x": 0, "y": 0},
            {"residue": "Gal", "x": 1, "y": 0},
            {"residue": "GalNAc", "x": 2, "y": 0},
        ],
        "edges": [(0, 1), (1, 2)],
    },
    "Core2": {
        "nodes": [
            {"residue": "Gal", "x": 0, "y": 0},
            {"residue": "GalNAc", "x": 1, "y": 0},
            {"residue": "GlcNAc", "x": 1, "y": -1},
        ],
        "edges": [(0, 1), (2, 1)],
    },
    "Core2_Fuc": {
        "nodes": [
            {"residue": "Gal", "x": 0, "y": 0},
            {"residue": "GalNAc", "x": 1, "y": 0},
            {"residue": "GlcNAc", "x": 1, "y": -1},
            {"residue": "Fuc", "x": 1, "y": 1},
        ],
        "edges": [(0, 1), (2, 1), (3, 1)],
    },
    "Core3": {
        "nodes": [
            {"residue": "GlcNAc", "x": 0, "y": 0},
            {"residue": "GalNAc", "x": 1, "y": 0},
        ],
        "edges": [(0, 1)],
    },
}


def _require_core_type(core_type: str) -> None:
    if core_type not in CORE_TYPE_COMPOSITION:
        raise ValueError(f"unknown core_type: {core_type}")


def glycoform_mass(core_type: str) -> float:
    """Return monoisotopic mass in Da. Unknown core_type -> ValueError."""
    composition = glycoform_composition(core_type)
    return sum(
        MONOISOTOPIC_RESIDUE_MASS[residue] * count
        for residue, count in composition.items()
    )


def glycoform_composition(core_type: str) -> dict:
    """Return a {monosaccharide: count} dict. Unknown -> ValueError."""
    _require_core_type(core_type)
    return dict(CORE_TYPE_COMPOSITION[core_type])


def glycoform_label(core_type: str) -> str:
    """Return a human-friendly label for UI tables."""
    _require_core_type(core_type)
    return CORE_TYPE_LABEL.get(core_type, core_type)


def _render_shape(residue: str, cx: float, cy: float, size: float) -> str:
    spec = SNFG_SPEC[residue]
    half = size / 2.0
    stroke = ' stroke="#000000" stroke-width="1"'
    if spec["shape"] == "square":
        x = cx - half
        y = cy - half
        return (
            f'<rect x="{x:.1f}" y="{y:.1f}" width="{size:.1f}" height="{size:.1f}"'
            f' fill="{spec["fill"]}"{stroke}/>'
        )
    if spec["shape"] == "circle":
        return (
            f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="{half:.1f}" fill="{spec["fill"]}"'
            f"{stroke}/>"
        )
    if spec["shape"] == "diamond":
        points = [
            (cx, cy - half),
            (cx + half, cy),
            (cx, cy + half),
            (cx - half, cy),
        ]
        points_attr = " ".join(f"{x:.1f},{y:.1f}" for x, y in points)
        return f'<polygon points="{points_attr}" fill="{spec["fill"]}"{stroke}/>'
    if spec["shape"] == "triangle":
        points = [
            (cx, cy - half),
            (cx + half, cy + half),
            (cx - half, cy + half),
        ]
        points_attr = " ".join(f"{x:.1f},{y:.1f}" for x, y in points)
        return f'<polygon points="{points_attr}" fill="{spec["fill"]}"{stroke}/>'
    raise ValueError(f"unsupported SNFG shape for {residue}: {spec['shape']}")


def render_glycoform_svg(
    core_type: str,
    symbol_px: int = 22,
    spacing_px: int = 8,
) -> str:
    """
    Render an SNFG-style SVG string for the glycoform.

    Layout: horizontal, reducing end on the LEFT (attached to protein).
    Glycoforms with branches (Core2 = Core1 + GlcNAc branch) render the
    branch above the main chain connected by a short vertical line.

    Returns a self-contained <svg>...</svg> string with no external deps
    (no CSS, no JS). Inline shapes only. Viewport sized to fit the shapes.
    """
    _require_core_type(core_type)
    if symbol_px <= 0:
        raise ValueError("symbol_px must be > 0")
    if spacing_px < 0:
        raise ValueError("spacing_px must be >= 0")

    layout = _LAYOUTS[core_type]
    step = float(symbol_px + spacing_px)
    padding = float(max(4, symbol_px // 3))
    y_step = float(symbol_px + max(4, spacing_px // 2))
    half = symbol_px / 2.0

    placed_nodes = []
    for node in layout["nodes"]:
        placed_nodes.append(
            {
                "residue": node["residue"],
                "cx": padding + node["x"] * step + half,
                "cy": padding + (node["y"] + 1) * y_step + half,
            }
        )

    width = max(node["cx"] for node in placed_nodes) + half + padding
    height = max(node["cy"] for node in placed_nodes) + half + padding

    parts = [
        (
            f'<svg xmlns="http://www.w3.org/2000/svg" width="{width:.0f}" '
            f'height="{height:.0f}" viewBox="0 0 {width:.1f} {height:.1f}" '
            'role="img" aria-label="glycan symbol">'
        )
    ]

    for start_idx, end_idx in layout["edges"]:
        start = placed_nodes[start_idx]
        end = placed_nodes[end_idx]
        parts.append(
            f'<line x1="{start["cx"]:.1f}" y1="{start["cy"]:.1f}" '
            f'x2="{end["cx"]:.1f}" y2="{end["cy"]:.1f}" '
            'stroke="#000000" stroke-width="1"/>'
        )

    for node in placed_nodes:
        parts.append(_render_shape(node["residue"], node["cx"], node["cy"], float(symbol_px)))

    parts.append("</svg>")
    return "".join(parts)
