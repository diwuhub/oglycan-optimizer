# oglycan-optimizer

> Physics-informed optimizer for EThcD acquisition parameters in O-glycopeptide mass spectrometry. Given a glycoprotein's site catalog, it tunes 6 instrument and workflow parameters to maximize site-localized glycan identification confidence. Pure Python stdlib, zero dependencies.

![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg) ![Python 3.9+](https://img.shields.io/badge/Python-3.9%2B-green.svg)

## Key results

| Metric | Value |
|--------|-------|
| Parameters optimized | 6 (SA%, MS2 resolution, HCD NCE, enzyme prep, LC gradient, search tolerance) |
| Sub-models | 6 literature-grounded Gaussian models composed into a single score |
| Reference catalogs | 5 glycoproteins (Etanercept, EPO, CTLA4-Ig, IgA1 hinge, Atacicept) |
| Dependencies | Zero (pure Python stdlib) |
| Uncertainty quantification | Bootstrap 95% CIs on composite score |
| Pilot-data loop | Ingest MSFragger-Glyco TSV, recalibrate site difficulties from observed localization |

## Problem

Every open-source glycoproteomics tool (GlycReSoft, Byonic, Skyline, pGlyco3, MSFragger-Glyco) answers: *given raw MS data, identify the glycopeptides.* None answer the question that comes first: *what acquisition parameters should I use to collect the raw data?*

For O-glycopeptide EThcD methods, the parameter space is unintuitive:

- ETD and HCD have **opposite** optimal SA% — ETD peaks at ~18%, HCD at ~35%
- Higher MS2 resolution improves mass accuracy but destroys scan rate (duty-cycle penalty)
- OpeRATOR enzyme requires SialEXO pretreatment to handle sialylated sites (+17% yield if included)
- Gradient length saturates above ~150 minutes

This tool composes six literature-grounded Gaussian sub-models into a single composite score, then lets you search parameter space against your specific glycoprotein's site catalog.

## Quick start

```bash
# No pip dependencies — just Python 3.9+
python -m oglycan tune sites/etanercept.json
```

Output: per-sub-model scores + composite score + sensitivity breakdown.

## CLI subcommands

| Subcommand | Description |
|------------|-------------|
| `tune` | Evaluate a site catalog and print the text report |
| `uncertainty` | Estimate a composite-score 95% confidence interval |
| `sensitivity` | Report signed local derivatives around the current method |
| `recommend` | Propose constrained next-run parameter moves |
| `ingest-pilot` | Normalize pilot evidence from MSFragger-Glyco TSV or canonical JSON |
| `recalibrate` | Update site difficulties from observed pilot localization |
| `suggest-catalog` | Build a starter site catalog from sequence candidates |
| `method-card` | Write a markdown method summary for manual transcription |
| `serve` | Launch the Streamlit UI (`pip install .[web]`) |

## Architecture

Six focused sub-models compose into one score:

```
src/oglycan/
├── models/
│   ├── etd_efficiency.py      — SA% ~18%, charge state coverage, isolation width
│   ├── hcd_fragments.py       — NCE ~30, glycan DB coverage, core-type specificity
│   ├── ms_acquisition.py      — Resolution vs duty cycle, injection time, AGC
│   ├── lc_separation.py       — Gradient slope, column class, saturation
│   ├── operator_bonus.py      — OpeRATOR + SialEXO prerequisite chain
│   └── search_quality.py      — MS1/MS2 mass tolerance, FDR threshold
├── site_scoring.py             — Per-site difficulty x glycoform x coverage
└── core.py                     — Composition: sub-models → site scores → composite
```

Each sub-model is <80 LOC, unit-tested independently, and its Gaussian optimum cites specific literature (see [docs/SCIENCE.md](docs/SCIENCE.md)).

## What this tool does NOT do

- **Does not process raw MS data.** GlycReSoft, Byonic, and MSFragger-Glyco solve that problem.
- **Does not identify specific glycan structures.** GlycoWorkbench and similar tools handle annotation.
- **Does not optimize sample prep beyond enzyme choice.** Sample-prep DOE is a separate concern.
- **Does not simulate full MS2 spectra.** Isotope modeling and fragmentation pathway detail are v2.0 goals.

See [LIMITATIONS.md](LIMITATIONS.md) for the full scope boundary.

## Who is this for

- **Mass spectrometry method developers** in biotech/pharma doing CMC-relevant glycoprotein characterization
- **Fc-fusion protein programs** (etanercept, abatacept, alefacept-class) where O-glycan site heterogeneity drives PK/PD variability
- **Bioprocess teams** studying how manufacturing changes affect O-glycoform distribution
- **Glycoproteomics methods groups** formalizing design-of-experiments for acquisition workflows

## Citation

```bibtex
@software{wu2026oglycanoptimizer,
  author = {Wu, Di},
  title  = {oglycan-optimizer: Physics-Informed EThcD Parameter Optimization for O-Glycoproteomics},
  year   = {2026},
  url    = {https://github.com/diwuhub/oglycan-optimizer}
}
```

## References

- Riley et al. 2020 Anal. Chem. — EThcD supplemental activation optimization for glycopeptides
- Malaker et al. 2022 Nat. Methods — OpeRATOR enzyme for O-glycoproteomics
- Saba et al. 2012 J. Proteome Res. — MS acquisition parameter effects on glycopeptide identification
- Struwe et al. 2017 Glycobiology — O-glycan site localization challenges

## License

MIT. See [LICENSE](LICENSE).
