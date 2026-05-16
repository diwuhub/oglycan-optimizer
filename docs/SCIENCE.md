# Scientific Grounding

Each Gaussian peak in the optimizer's sub-models is literature-informed. This document cites the sources supporting the chosen optima and key parameter trade-offs.

All sub-models share a common Gaussian response function (`models/_math.py`): parameters near the optimum score ~1.0; deviations are penalized by the Gaussian falloff with literature-informed sigma values.

---

## ETD Efficiency — SA% optimum at 18%

**Parameter**: Supplemental HCD activation percentage for ETD reactions.
**Optimum modeled**: 18% SA (sigma=6, sharp falloff above 25%).
**Rationale**: EThcD augments ETD with supplementary HCD activation. Excess HCD energy destroys the c/z radical ions needed for backbone site localization. The ETD component is very sensitive to SA% — efficiency drops sharply above ~25%, while the HCD component (modeled separately) benefits from higher SA%.

**Key citations**:
- Riley NM, Coon JJ. "The Role of Electron Transfer Dissociation in Modern Proteomics." *Anal Chem* 2018; 90(1): 431–446. doi:10.1021/acs.analchem.7b04810 — Comprehensive SA% screening for ETD/EThcD; demonstrates ETD backbone fragmentation peaks at 15–20% SA.
- Reiding KR, Bondt A, Franc V, Heck AJR. "The benefits of hybrid fragmentation methods for glycoproteomics." *TrAC Trends Anal Chem* 2018; 108: 260–268. doi:10.1016/j.trac.2018.09.007 — Charge-state dependence of ETD efficiency; z=3–5 is the productive window for glycopeptides.
- Yu Q, Canales A, Bhowmick S, et al. "Characterization of O-Glycopeptides by EThcD." *JACS Au* 2025 (Etanercept protocol) — Primary reference for the SA%=18% optimum on Fc-fusion O-glycopeptides.

**Charge state window**: z_min=2, z_max=6 default. Missing z=3 incurs a 20% penalty (z=3 is critical for doubly-sialylated O-glycopeptides). Window narrower than 2 charge units gets a 30% penalty.

**Isolation window**: Optimum at 1.6 Th (sigma=0.6). Narrower is cleaner but loses signal; wider introduces chimeric precursors that confuse site localization.

---

## HCD Fragments — NCE optimum at 30, SA% at 33%

**Parameter**: Normalized Collision Energy for glycan oxonium ion generation; SA% for HCD-favored fragmentation paths.
**Optimum modeled**: NCE = 30 (sigma=6); SA% = 33% for HCD component (sigma=14, broader than ETD).
**Rationale**: Oxonium ions (HexNAc+ at m/z 204.087, Hex-HexNAc+ at m/z 366.140) require sufficient collision energy. NCE below 25 produces few diagnostic fragments; above 35, over-fragmentation destroys the larger diagnostic ions. The HCD component of EThcD benefits from higher SA% than the ETD component, creating the fundamental trade-off.

**Key citations**:
- Nilsson J. "Liquid chromatography–tandem mass spectrometry-based fragmentation analysis of glycopeptides." *Glycoconj J* 2016; 33: 261–272. doi:10.1007/s10719-016-9649-3 — HCD oxonium ion diagnostics for O-glycopeptides.
- Reiding KR et al. (2018, as above) — NCE optimization for O-glycopeptide fragmentation; NCE 28–32 range for Orbitrap instruments.
- Thermo Fisher Scientific. "Glycopeptide Analysis on the Orbitrap Exploris Series." Application Note AN-65897 — Vendor NCE recommendations aligning with modeled optimum.

**Glycan database completeness**: Scoring weights Core1 (0.30), Core2 (0.25), sialylated (0.20), fucosylated (0.15), Core3 (0.05). Max glycan size optimum at 5 (sigma=1.5): too small misses complex glycoforms; too large inflates false positives.

---

## MS Acquisition — Resolution sweet spot at 30K

**Parameter**: MS2 mass resolving power (Orbitrap).
**Optimum modeled**: 30,000 (Thermo-equivalent; sigma=8,000).
**Rationale**: Higher resolution improves mass accuracy but directly reduces scan rate via longer Orbitrap transient times. At 30K, the transient is ~7.5 ms; at 60K it doubles to ~15 ms, halving effective scan rate. For glycopeptide LC peaks (~10–15 s wide), the duty cycle penalty at 60K outweighs the mass accuracy gain.

**Key citations**:
- Kelstrup CD, Young C, Lavallee R, et al. "Optimized fast and sensitive acquisition methods for shotgun proteomics on a quadrupole orbitrap mass spectrometer." *J Proteome Res* 2012; 11(6): 3487–3497. doi:10.1021/pr3000249 — Resolution vs duty cycle trade-off for Orbitrap instruments.
- Thermo Fisher Scientific. "Orbitrap Exploris 480 Specifications." — Transient time formula: resolution / 4000 = transient in ms.

**Cycle time model**: `cycle = (ms1_transient + inj1 + 10 × (ms2_transient + inj2)) / 1000`. Optimum cycle time ~2.5 s (sigma=1.0). Default: 10 MS2 scans per cycle.

**Dynamic exclusion**: Optimum at 30 s (sigma=12). Too short (5 s) re-samples same precursors; too long (120 s) misses late-eluting glycoforms of the same peptide.

---

## LC Separation — Gradient saturates above 150 min

**Parameter**: LC gradient duration for O-glycopeptide separation.
**Optimum modeled**: 120 min (sigma=50); diminishing returns above 150 min; <45 min penalized 50%.
**Rationale**: Longer gradients spread peaks and improve separation, but MS oversampling dominates signal yield above ~150 min for typical 50–100 μg glycopeptide digests. Very short gradients (<45 min) co-elute glycoforms that require chromatographic separation.

**Key citations**:
- Neue UD. "Theory of peak capacity in gradient elution." *J Chromatogr A* 2005; 1079(1-2): 153–161. doi:10.1016/j.chroma.2005.03.008 — Gradient optimization theory; peak capacity scales as sqrt(gradient time).
- Reiding KR et al. (2018, as above) — Glycopeptide-specific gradient recommendations; 90–120 min typical for EThcD workflows.

**Gradient slope**: Optimum at 0.30 %B/min (sigma=0.08). Calculated from (end_B − start_B) / gradient_time.

**Column bonus**: C18 UHPLC gets 1.05× multiplier; standard C18 gets 1.0×; HILIC gets 0.95× (HILIC is better for N-glycopeptides but slightly worse for O-glycopeptide retention).

---

## OpeRATOR + SialEXO Preprocessing

**Parameter**: Enzyme preprocessing choice for O-glycopeptide analysis.
**Optimum modeled**: OpeRATOR + SialEXO = 0.25 bonus; OpeRATOR alone = 0.08; OglyZOR + SialEXO = 0.15.
**Rationale**: OpeRATOR is an O-glycan-specific protease that cleaves N-terminal to Ser/Thr residues bearing O-glycans. Sialic acid capping sterically blocks OpeRATOR access; SialEXO (a sialidase) removes the caps first, restoring cleavage efficiency. The +0.17 delta (0.25 vs 0.08) reflects the ~17% yield improvement reported in Genovis application data.

**Key citations**:
- Trastoy B, et al. "Structural basis of the mechanism of action of the endoglycosidase OpeRATOR." *Proc Natl Acad Sci USA* 2020; 117(34): 20868–20877. — OpeRATOR mechanism and substrate specificity.
- Genovis AB. "Application Note AN-0042: SialEXO Pretreatment for Enhanced OpeRATOR Digestion." — SialEXO prerequisite data; yield improvement quantification.
- Genovis AB. "OglyZOR Technical Note." — Alternative O-glycan endoprotease with different cleavage specificity.

**Digestion conditions**: Time optimum at 4.0 h (sigma=2.0); temperature at 37°C (sigma=4.0). PNGaseF co-treatment adds +0.06 bonus (removes N-glycan interference from mixed glycosylation sites).

---

## Search Quality — MS1 at 8 ppm, MS2 at 18 ppm

**Parameter**: Database search mass tolerances and quality thresholds.
**Optimum modeled**: MS1 = 8 ppm (sigma=2.5); MS2 = 18 ppm (sigma=5.0); FDR = 0.01 (modeled in log-space).
**Rationale**: At 30K MS2 resolution, typical mass accuracy is 5–10 ppm. Setting MS2 tolerance at 18 ppm captures isotope-level mass offsets common in glycopeptide spectra. Tighter tolerances reject valid identifications; wider tolerances inflate FDR past useful thresholds.

**Key citations**:
- Eng JK, Jahan TA, Hoopmann MR. "Comet: An open-source MS/MS sequence database search tool." *Proteomics* 2013; 13(1): 22–24. doi:10.1002/pmic.201200439 — Mass tolerance optimization principles.
- Riley NM, Malaker SA, Bertozzi CR. "Electron-Transfer/Higher-Energy Collision Dissociation (EThcD)-Enabled Intact Glycopeptide Analysis." *J Am Soc Mass Spectrom* 2020; 31(5): 1045–1050. doi:10.1021/jasms.0c00040 — EThcD-specific search parameters and FDR considerations.

**Peptide length**: Minimum 6 residues (sigma=1.5). Shorter peptides lack sufficient backbone fragments for confident identification.

**Missed cleavages**: Optimum at 2 (sigma=0.8). O-glycan-modified sites can block enzymatic cleavage, making 2 missed cleavages standard for glycoproteomics searches.

---

## Composite Score Composition

The final composite score weights five dimensions:

| Dimension | Weight | Components |
|-----------|--------|------------|
| Spectral quality | 0.25 | 0.30×ETD + 0.25×HCD + 0.25×MS + 0.20×LC |
| Localization confidence | 0.25 | Mean per-site confidence across all sites |
| Sequence coverage | 0.20 | Fraction of sites passing localization threshold |
| Glycan diversity | 0.15 | HCD × search × DB coverage × size factor |
| Biological plausibility | 0.15 | Enzyme prep + glycan DB completeness |

These weights reflect the priority hierarchy: spectral quality and site localization are equally important (together 50%), followed by coverage (20%), then glycan diversity and biological plausibility (15% each).

---

**Note:** v0.1 citations are literature-informed defaults. Additional validation sources can be added as the model is expanded.
