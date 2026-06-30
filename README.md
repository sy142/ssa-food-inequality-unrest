# Does income inequality moderate the effect of food price inflation on social unrest in Sub-Saharan Africa?

Replication materials for the manuscript *"Does income inequality moderate the effect of food price inflation on social unrest in Sub-Saharan Africa?"*

## Overview

This repository contains the data and code to reproduce the analysis. The study examines whether income inequality moderates the relationship between food price inflation and social unrest across Sub-Saharan African countries (2000–2025), using two-way fixed-effects count models (negative binomial and Poisson pseudo-maximum likelihood).

The central finding is methodological as well as substantive: an apparent positive moderation under conventional model-based inference does not survive country-clustered or pairs-cluster bootstrap inference. The repository allows full reproduction of the estimates, the inference comparison, the diagnostics, and the figures.

## Repository structure

```
.
├── README.md
├── LICENSE
├── CITATION.cff
├── data/
│   ├── ssa_food_unrest_panel.xlsx   processed analytic panel
│   └── DATA_SOURCES.md              variable sources and access notes
├── code/
│   ├── 01_build_panel.py            assembles the panel from source series
│   ├── 02_analysis.R                estimation, inference, diagnostics
│   └── 03_figures.py                figures
└── figures/
    ├── fig1_unrest_overview.png
    ├── fig2_marginal_effect.png
    └── fig3_ssa_inequality_arrows.png
```

## Requirements

R (4.5.2 or later) with: `plm`, `readxl`, `dplyr`, `car`, `glmmTMB`, `splines`, `fixest`, `sandwich`, `clubSandwich`, `lmtest`, `marginaleffects`, `modelsummary`, `MASS`, `boot`, `sf`, `rnaturalearth`, `rnaturalearthdata`, `ggplot2`, `scales`, `ggrepel`.

Python (3.11 or later) with: `pandas`, `numpy`, `country_converter`, `wbgapi`, `faostat`, `requests`, `openpyxl`, `matplotlib`.

## Usage

Run scripts from the repository root so that the relative paths (`data/`, `figures/`) resolve correctly.

1. `code/01_build_panel.py` assembles the country-year panel from the source series. Raw source data must be downloaded separately (see `data/DATA_SOURCES.md`); the script writes the processed panel as an Excel file.
2. `code/02_analysis.R` reproduces the estimation, the model-based versus cluster-robust and bootstrap inference, the panel diagnostics (cross-sectional dependence, CIPS unit-root, Dumitrescu–Hurlin causality), the marginal effects, the leave-one-country-out jackknife, the Southern-Africa exclusion, and the nutritional-vulnerability tests. It also writes intermediate files used by the figures.
3. `code/03_figures.py` reproduces Figure 1 (descriptive overview) and Figure 2 (marginal effect with model-based and country-clustered confidence bands). Figure 3 (the map of country-specific marginal effects) is produced inside `02_analysis.R` using `ggplot2`, `sf` and `rnaturalearth`.

The processed panel (`data/ssa_food_unrest_panel.xlsx`) is provided so that the analysis and figures can be run directly without re-running the build step.

## Data availability

The processed analytic panel is provided here. Raw source data are not redistributed; they are publicly available from the original providers listed in `data/DATA_SOURCES.md`. ACLED data are subject to the ACLED terms of use and must be obtained directly from ACLED.

## Citation

If you use these materials, please cite the manuscript and this repository (see `CITATION.cff`).

## Author

Salim Yılmaz (ORCID: [0000-0003-2405-5084](https://orcid.org/0000-0003-2405-5084))

## License

Code is released under the MIT License (see `LICENSE`). Data files are subject to the licenses of their original providers.
