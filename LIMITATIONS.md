# Limitations

Honest disclosure of v0.1 scope boundaries.

## What the v0.1 MVP does NOT model

### No isotope peak overlap simulation

Real glycopeptides produce overlapping isotope clusters. At any given MS2 resolution, nearby peaks merge, reducing site localization confidence. The current engine scores resolution as a single Gaussian; it does not compute actual peak-recognition probability given specific site mass-spacing.

**Consequence:** The optimizer may over-recommend low MS2 resolution for glycoproteins with closely-spaced O-glycosites (e.g., Etanercept's T237/T239 stalk cluster). Known gap; v2.0 target.

### No fragmentation pathway detail

Oxonium ions, HexNAc-loss patterns, Gal-loss patterns, and Core1-vs-Core2-specific fragments each have distinct formation kinetics. The current engine reduces fragmentation to one Gaussian (SA% → oxonium intensity). It does not distinguish which fragment types are most diagnostic for site localization vs glycan structure determination.

**Consequence:** For glycoproteins dominated by a single core type, the scoring may miss opportunities to optimize for that core's characteristic fragments. v2.1 target.

### No kinetic model for enzyme digestion

OpeRATOR digestion is scored on time + temperature + SialEXO prerequisite. A full Michaelis-Menten kinetic model with substrate accessibility would better predict yield for low-abundance sites.

**Consequence:** For sites with unusually hindered substrate access, the tool may under-estimate the benefit of prolonged digestion. v2.2 target.

### No cell-culture sensitivity modeling

Difficulty scores for each site are assigned from published observations. Real glycoform distributions shift with bioprocess parameters (pH, dissolved oxygen, feed strategy). The tool does not predict these shifts.

**Consequence:** Method optimized for bioreactor run A may be suboptimal for run B with different bioprocess conditions. v2.2 target.

### Gaussian-only parameter surfaces

All six sub-models assume Gaussian optima. Real MS response functions can be multimodal (e.g., high vs low charge states have different optimal conditions) or flat-topped (saturation plateaus). Gaussian fits are pedagogically clear and literature-grounded for central behavior, but miss edge-case topology.

**Consequence:** Near optima, predictions are reliable; far from optima, predictions may mislead. Users should inspect sub-model plots (future feature) to confirm working point is in the well-fit region.

### Etanercept-specific site catalog in v0.1

The reference implementation is Etanercept's 13 O-glycosites. Additional glycoproteins (EPO, CTLA4-Ig, IFN-beta-Fc) are v1.1 deliverables. Users can curate their own site catalogs now following `sites/schema.md`, but should validate their curation against published literature before trusting optimizer output.

## What Gaussian parameters are based on

The optimal values used in v0.1 come from:

- **JACS Au 2025 Etanercept EThcD protocol** — primary reference for SA% optimum, enzyme workflow
- **Riley et al., Anal Chem 2020** — EThcD supplemental activation screen
- **Orbitrap instrument specification documents** — duty cycle and transient time formulas
- **OpeRATOR / SialEXO product literature** — enzyme kinetics ranges

These are **informed defaults**, not experimentally validated constants for arbitrary glycoproteins. Every citation is explicit in [docs/SCIENCE.md](docs/SCIENCE.md).

## What the tool assumes about the instrument platform

v0.1 math is calibrated for Orbitrap-class instruments (Fusion / Exploris / Astral family). It does NOT automatically adjust for:

- Timstof-class (ion-mobility) instruments
- Bruker instruments (different transient-time scaling)
- Vendor-specific AGC behaviors

Porting to other platforms requires adjusting constants in `models/ms_acquisition.py`. Multi-platform support is outside the current package scope.

## Composite score interpretability

The composite score (0-1 scale) is a weighted combination of six sub-model outputs. It is useful for ranking parameter sets against each other (A is better than B for this glycoprotein). It is NOT an absolute predictor of identification rate — a score of 0.95 does not mean 95% of sites will be successfully localized in real data.

For absolute performance prediction, you still need to run acquisition on a representative sample.

## No warranty

Per the MIT license, this software is provided "as is" without warranty. Scientific conclusions drawn from optimized methods must be validated experimentally.

## Reporting issues

Found a case where the optimizer contradicts your experimental results? Open a GitHub issue with:

1. Your glycoprotein and site catalog
2. The parameter set the optimizer recommended
3. What worked better experimentally and why you think so

This evidence helps calibrate parameter confidence.

## v1.1 empirical findings

### Composite CI is wide

Monte Carlo with literature-informed `mu` uncertainties (`sigma / 3`) gives the default Etanercept method a 95% CI of roughly `[0.68, 0.94]`, for a width of about `0.26`.

The point estimate sits at the Gaussian peaks, so perturbations can only lower the composite. In practice the sampled distribution is one-sided below the point estimate rather than symmetric around it.

### `recommend_next_method` is coordinate descent, not joint optimization

The recommender does a single deterministic pass over a hardcoded 6-parameter search space.

That keeps it fast and reproducible, but it will miss improvements that require simultaneous moves in multiple parameters. Joint search is outside the current package scope.

### Cross-protein generalization check is preliminary

All 5 reference catalogs score at least `0.93` under the defaults, so the current floor-test (`>= 0.50`) does not discriminate between easy and brittle settings.

A stronger robustness check would stress the parameters, then assert that relative ranking stays consistent across proteins.

### Method export is a markdown card, not Xcalibur XML

Thermo method files are version-specific and proprietary.

Emitting synthetic XML risks generating files that load incorrectly in the method editor, so the exporter writes a markdown card for manual transcription instead.

### No real predictor is shipped

`suggest-catalog --predictor scan_st` returns every Ser/Thr as a candidate with `p_glycosite=null`.

NetOGlyc and Stack-OglyPred-PLM adapters currently raise `NotImplementedError` with guidance strings; real integration is future opt-in work that would add external dependencies.

### Only MSFragger-Glyco `psm.tsv` has a real pilot parser

Byonic and O-Pair import paths are still stubs.

### No wet-lab validation yet

The composite score is a ranking metric for comparing candidate methods.

Absolute performance against real acquisitions has not yet been measured on experimental runs using tool-recommended parameters.
