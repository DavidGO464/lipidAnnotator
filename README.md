# lipidAnnotator <img src="https://img.shields.io/badge/R-package-blue" align="right"/>

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![LIPID MAPS](https://img.shields.io/badge/nomenclature-LIPID%20MAPS%202020-green)](https://www.lipidmaps.org)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19325471.svg)](https://doi.org/10.5281/zenodo.19325471)

> Universal lipid name parser and annotator for R, with canonical LIPID MAPS
> classification and shorthand notation per Liebisch et al. 2020.

---

## Overview

Lipid nomenclature in mass spectrometry-based lipidomics is inconsistent across
software platforms (LipidSearch, MS-DIAL, LipidView). **lipidAnnotator** provides
a single function — `annotate_lipid()` — that parses any lipid name format and
returns a structured data frame with LIPID MAPS canonical classification,
chain-level metadata, and optional shorthand notation.

### Supported input formats

| Input | Class |
|---|---|
| `PC 16:0/18:1` | Diacyl GPL |
| `PC O-18:1/20:4` | Ether GPL |
| `PE P-18:0/22:6` | Plasmalogen |
| `PE-Cer 20:1/16:0` | Sphingo-GPL hybrid |
| `Cer(d18:0/C22:0)` | Ceramide (legacy C-prefix) |
| `TG(16:0/18:1/18:1)` | Triacylglycerol |
| `HexCer 18:1/16:0` | Hexosylceramide |
| `Lyso PE 18:1(d7)` | Lyso-GPL + ISTD |
| `CL 18:1/18:1/18:1/18:1` | Cardiolipin (4 chains) |
| `FAHFA 18:1/9-O-16:0` | Fatty acid estolide |
| `15-MHDA` | Trivial name → FA 17:0 |
| `15-HETE` | Eicosanoid mediator |
| `Resolvin D1` | Specialized pro-resolving mediator |
| `WE 16:0/18:1` | Wax ester |
| `CoA 16:0` | Acyl-CoA |

---

## Installation
```r
# Install from GitHub
devtools::install_github("DavidGO464/lipidAnnotator")
```

---

## Usage
```r
library(lipidAnnotator)

lipids <- c(
  "PC O-18:1/20:4", "PE P-18:0/22:6",
  "PE-Cer 20:1/16:0", "HexCer 18:1/16:0",
  "TG(16:0/18:1/18:1)", "Cer d18:1/16:0",
  "Lyso PE 18:1(d7)", "CL 18:1/18:1/18:1/18:1",
  "15-HETE", "Resolvin D1", "WE 16:0/18:1", "CoA 16:0"
)

# Compact — equivalent to lipidR output (default)
annotate_lipid(lipids)

# Standard — adds LIPID MAPS classification and structural flags
annotate_lipid(lipids, detail = "standard")

# Full — all columns, for join with LMSD or reporting
annotate_lipid(lipids, detail = "full", shorthand = TRUE)
```

### Output levels

| Column | compact | standard | full |
|---|:---:|:---:|:---:|
| `Class`, chains, totals | ✓ | ✓ | ✓ |
| `lm_category`, `lm_class_id` | | ✓ | ✓ |
| `is_ether`, `is_plasmalogen`, `is_istd` | | ✓ | ✓ |
| `lm_class_name`, `headgroup`, `clean_name` | | | ✓ |
| `hydroxyl_base`, `ether_link`, `class_stub` | | | ✓ |
| `shorthand_lm` (Liebisch 2020) | optional | optional | optional |

### Key features

- **Ether / plasmalogen detection** — `PC O-18:1/20:4` → `is_ether = TRUE`; `PE P-18:0/22:6` → `is_plasmalogen = TRUE`
- **Sphingoid prefix** — `d18:1` correctly parsed as dihydroxy base, not flagged as ISTD
- **ISTD detection** — `18:1(d7)` → `is_istd = TRUE`; `d18:1` → `is_istd = FALSE`
- **Shorthand notation** — canonical output per Liebisch et al. 2020: `HexCer 18:1/16:0` → `HexCer(d18:1/16:0)`
- **LIPID MAPS ontology** — category (FA/GL/GP/SP/ST), class ID (e.g. `GP0101`), class name, headgroup
- **Lipid mediator parsing** — trivial names with positional notation: `15-HETE`, `12(13)-EpOME`, `Resolvin D1`, `PGE2`, `LTB4`

---

## Supported lipid classes

Nomenclature follows [LIPID MAPS](https://www.lipidmaps.org/) conventions
(Liebisch et al. 2020; Conroy et al. 2024). Ether linkages (O-/P-),
plasmalogens, deuterium-labeled ISTDs, and sphingoid base prefixes (d/m/t)
are detected automatically across all supported classes.

| Category | Classes |
|---|---|
| Fatty Acyls (FA) | FA, FAHFA, FAL, FOH, CAR, NAE, WE, CoA · Eicosanoids: HETE, HEPE, EET, DHET, OxoETE, HpETE, EpETE, HHTrE · Prostanoids: PGE2, PGD2, PGF2a, PGI2, PGB2, TXB1/2/3 · Leukotrienes: LTB4, LTC4, LTD4 · Octadecanoids: HODE, HOTrE, EpOME · SPMs: Resolvins, Maresins, Protectins, Lipoxins, HDoHE |
| Glycerolipids (GL) | MG, DG, TG, MGDG, DGDG, SQDG, SQMG |
| Glycerophospholipids (GP) | PA/LPA, PC/LPC, PE/LPE, PG/LPG, PI/LPI, PS/LPS, PIP/PIP2/PIP3, CL, MLCL, DLCL, BMP, LBPA, PEtOH, PThr, plasmenylPC/PE |
| Sphingolipids (SP) | Cer, CerP, SM, LSM, HexCer, LHexCer, Hex2Cer, GlcCer, GalCer, LacCer, SHexCer, GM1/2/3, GD1/GD3, GT3, Sph, S1P, Sa, So, Sa1P, So1P, PE-Cer, PI-Cer, PA-Cer, IPC, deoxyCer |
| Sterol Lipids (ST) | Chol, CE, FC, BA |
| Prenol Lipids (PR) | CoQ |

---

## Integration with lipidomics pipelines
```r
# Annotate features from a limma result
lipid_annot <- annotate_lipid(rownames(limma_results), detail = "standard")

results_annotated <- merge(
  limma_results,
  lipid_annot[, c("Molecule", "Class", "lm_category",
                  "lm_class_id", "is_ether", "total_cl", "total_cs")],
  by.x = "row.names", by.y = "Molecule"
)
```

---

## References

- Liebisch G et al. Update on LIPID MAPS classification, nomenclature, and shorthand notation for MS-derived lipid structures. *J Lipid Res.* 2020;61(12):1539–1555. [PMID: 33037133](https://pubmed.ncbi.nlm.nih.gov/33037133/)
- Conroy MJ et al. LIPID MAPS: update to databases and tools for the lipidomics community. *Nucleic Acids Res.* 2024;52(D1):D1677–D1682. [PMID: 37855672](https://pubmed.ncbi.nlm.nih.gov/37855672/)
- Fahy E et al. Update of the LIPID MAPS comprehensive classification system for lipids. *J Lipid Res.* 2009;50 Suppl:S9–S14. [PMID: 19098281](https://pubmed.ncbi.nlm.nih.gov/19098281/)

---

## License

MIT ©
