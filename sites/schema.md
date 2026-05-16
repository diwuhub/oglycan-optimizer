# Glycoprotein Site Catalog Schema

Each glycoprotein file in this directory describes its O-glycosylation site catalog for the optimizer.

## Schema

```json
{
  "glycoprotein": {
    "name": "string",
    "description": "string",
    "reference": "string (citation)",
    "uniprot": "string (optional UniProt accession)"
  },
  "sites": [
    {
      "pos": "integer (residue number)",
      "aa": "string (amino acid, S or T)",
      "domain": "string (e.g., 'stalk', 'CRD', 'Fc')",
      "difficulty": "float [0.4, 0.9] (intrinsic localization difficulty)",
      "core_types": ["Core1", "Core1_Sia", "Core2", "Core2_Fuc"],
      "immunogenicity_flags": ["sialylation_PK", "fucosylation_ADCC"]
    }
  ],
  "predicted_sites": [
    {
      "pos": "integer (residue number)",
      "aa": "string (amino acid, S or T)",
      "p_glycosite": "float [0.0, 1.0] or null",
      "source": "string (e.g., 'scan_st', 'netoglyc_4.0')"
    }
  ],
  "n_glycosites": [
    {
      "pos": "integer (residue number)",
      "domain": "string"
    }
  ],
  "localization_threshold": "float (0.0-1.0)"
}
```

## Fields (O-glycosites)

- **pos**: residue number (1-indexed, per canonical sequence)
- **aa**: amino acid type at this position (S or T, the only O-glycosylation substrates)
- **domain**: structural context (e.g., 'stalk', 'hinge', 'CH2', 'CRD', 'Fc')
- **difficulty**: 0.4 = well-resolved in standard methods; 0.9 = requires bespoke optimization. Derived from prior experimental evidence.
- **core_types**: observed core structures at this site (Core1, Core2, and their modified variants like Core1_Sia, Core2_Fuc)
- **immunogenicity_flags**: structural features with clinical/CMC risk (e.g., sialylation_PK, fucosylation_ADCC)

## Fields (N-glycosites)

- **pos**: residue number (part of NxS/T motif)
- **domain**: structural context

## Fields (Predicted Sites, optional)

- **predicted_sites**: review-only candidate O-glycosites. These are not scored by the optimizer until a user promotes experimentally supported observations into `sites`.

## Fields (Metadata)

- **localization_threshold**: confidence threshold for accepting localization assignments (typically 0.75)

## Adding a New Glycoprotein

1. Copy this schema
2. Curate sites from literature or internal MS experiments
3. Save as `{glycoprotein_slug}.json`
4. Run: `python -m oglycan tune sites/{glycoprotein_slug}.json`

## Available Catalogs (MVP)

- `etanercept.json` — TNFR2-Fc fusion, 13 O-sites (reference implementation)

## Additional Example Catalogs

- `epo.json` — Erythropoietin O-sites
- `ctla4_ig.json` — CTLA4-Ig fusion
- `ifn_beta_fc.json` — Interferon-beta-Fc
