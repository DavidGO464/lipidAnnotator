
# lipidAnnotator

An R package for parsing lipid names and annotating them with LIPID MAPS 
canonical classification and shorthand notation.

## Features

- Parses any lipid name format (diacyl, ether, lyso, sphingolipids, neutral lipids)
- LIPID MAPS classification (category, class ID, class name)
- Canonical shorthand notation per Liebisch et al. J Lipid Res 2020
- Three output detail levels: compact, standard, full
- Correct detection of ether/plasmalogen linkages (O-/P-)
- Sphingoid base prefix handling (d/m/t)
- ISTD detection without false positives on sphingoid prefixes

## Installation

```r
# Install from GitHub
devtools::install_github("DavidGO464/lipidAnnotator")
```

## Example

This is a basic example which shows you how to solve a common problem:

```r
library(lipidAnnotator)

lipids <- c("PC O-18:1/20:4", "PE-Cer 20:1/16:0", "HexCer 18:1/16:0",
            "TG(16:0/18:1/18:1)", "Cer d18:1/16:0", "Lyso PE 18:1(d7)")

# Compact output (default)
annotate_lipid(lipids)

# With LIPID MAPS classification
annotate_lipid(lipids, detail = "standard")

# Full output with canonical shorthand
annotate_lipid(lipids, detail = "full", shorthand = TRUE)
```


## References

- Liebisch G et al. J Lipid Res. 2020;61(12):1539-1555. PMID: 33037133
- Conroy MJ et al. Nucleic Acids Res. 2024;52(D1):D1677-D1682. PMID: 37855672
- Fahy E et al. J Lipid Res. 2009;50 Suppl:S9-14. PMID: 19098281

