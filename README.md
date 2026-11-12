# Spatial Data Science and Disease Mapping in R

Modernized implementations, statistical theory, and diagnostic workflows for areal spatial data science, spatial autocorrelation, and hierarchical Bayesian disease mapping in R.

---

## Notebooks

| Notebook | Focus | Primary Dataset | Key Concepts & Methods |
| :--- | :--- | :--- | :--- |
| [**`main.ipynb`**](file:///home/olande/Documents/GitHub/Spatial/main.ipynb) | Nationwide Cancer Mortality ESDA | US Contiguous Counties (`cancer_reg.csv`, $n=3,022$) | Deterministic FIPS harmonization, age-standardization orthogonality test, Gaussian moments vs. Poisson overdispersion ($\text{Var}/\text{E} \approx 1,553$), small-area funnel variance, bivariate screening under FWER control (Bonferroni & Holm), collinearity diagnostics (VIFs $\approx 5.0$, condition number $\kappa = 6.30$), Queen contiguity graph export (`us_counties.graph`), Global Moran's $I$ ($0.462$), OLS residual autocorrelation ($0.342, p < 10^{-200}$), Lagrange Multiplier diagnostics ($\text{Robust LM Error} \gg \text{Robust LM Lag}$), 6-lag correlogram, and LISA cluster maps under FDR ($q < 0.05$) and Bonferroni thresholds. |
| [**`Spatial-Neighborhood-Matrices.ipynb`**](file:///home/olande/Documents/GitHub/Spatial/Spatial-Neighborhood-Matrices.ipynb) | Spatial Adjacency & Topology | Columbus, OH (`columbus.gpkg`, $n=49$) | Queen vs. Rook contiguity, $k$-NN, distance bands ($\le 0.62$), higher-order topological lags (`nblag`), row-standardization ($W$), and inverse-distance gravity weights. |
| [**`Spatial-Autocorrelation.ipynb`**](file:///home/olande/Documents/GitHub/Spatial/Spatial-Autocorrelation.ipynb) | Global & Local Autocorrelation | Boston, MA (`boston_tracts.gpkg`, $n=506$) | Global Moran's $I$, Monte Carlo permutation distributions ($M=999$), Moran scatterplot regression equivalence, Local Moran ($I_i$), and Anselin LISA cluster maps (Hotspots vs. Coldspots). |
| [**`Bayesian-Spatial-Models.ipynb`**](file:///home/olande/Documents/GitHub/Spatial/Bayesian-Spatial-Models.ipynb) | Hierarchical Spatial Regression | Boston, MA (`boston_tracts.gpkg`, $n=506$) | Fast deterministic inference with INLA, classical BYM vs. modern scaled BYM2 models, spatial precision $\tau_s$, mixing parameter $\phi$, and non-linear dollar back-transformations (`inla.tmarginal`). |
| [**`Disease-Risk-Modeling.ipynb`**](file:///home/olande/Documents/GitHub/Spatial/Disease-Risk-Modeling.ipynb) | Small-Area Epidemiological Risk | Pennsylvania (`pennLC`, $n=67$ counties) | Demographic indirect standardization (16 strata), expected counts ($E_i$), raw Standardized Mortality Ratios (SMR), Poisson log-linear spatial models, Bayesian shrinkage, and exceedance probabilities ($P(RR > 1.20)$). |
| [**`Areal-Data-Issues.ipynb`**](file:///home/olande/Documents/GitHub/Spatial/Areal-Data-Issues.ipynb) | Methodological Hazards & Solutions | 256 simulated micro-units & regional strata | The Modifiable Areal Unit Problem (MAUP scale and zoning effects), the Ecological Fallacy (Robinson's paradox), spatial change of support (COSP), Bayesian melding, and downscaling. |

---

## Datasets

The `data/` directory contains clinical oncology and census demographic datasets alongside spatial topology files:

* **`data/cancer_reg.csv`**: County-level age-adjusted cancer mortality rates (`TARGET_deathRate`) combined with 33 demographic, economic, and clinical covariates compiled from SEER and the US Census American Community Survey (ACS).
* **`data/death.csv` & `data/incd.csv`**: CDC/NCI SEER mortality and incidence tables with validated 5-digit county Federal Information Processing Standards (FIPS) codes.
* **`data/us_counties.graph`**: Pre-computed topological adjacency graph covering 3,022 contiguous US counties in binary format, ready for direct loading in spatial modeling software via `inla.read.graph()`.

---

## Technical Standards

* **Modern R & Tidyverse**: All data transformations use modern Tidyverse idioms with the native pipe (`|>`). Base R legacy assignment and filtering idioms are strictly avoided.
* **Package Management**: Managed through `librarian::shelf(...)` (`sf`, `spdep`, `tigris`, `INLA`, `SpatialEpi`, `tidyverse`, `patchwork`, `viridis`).
* **Visualizations**: Multi-panel figures composed with `patchwork`. All plotting blocks specify dual-kernel canvas geometry (`ark.plot.*` for Positron and `repr.plot.*` for IRkernel/Jupyter) for crisp, scrollbar-free rendering.
* **Mathematical & Statistical Rigor**: Statistical workflows include formal hypothesis tests, multiple testing corrections (Bonferroni, Holm, False Discovery Rate), variance inflation factors, and residual spatial diagnostics.
