# Spatial Data Science and Disease Mapping in R

Modernized implementations, statistical theory, and diagnostic workflows for areal spatial data science, spatial autocorrelation, and hierarchical Bayesian disease mapping in R.

---

## Notebooks

| Notebook | Focus | Primary Dataset | Key Concepts & Methods |
| :--- | :--- | :--- | :--- |
| [**`Spatial-Neighborhood-Matrices.ipynb`**](file:///home/olande/Documents/GitHub/Spatial/Spatial-Neighborhood-Matrices.ipynb) | Spatial Adjacency & Topology | Columbus, OH (`columbus.gpkg`, $n=49$) | Queen vs. Rook contiguity, $k$-NN, distance bands ($\le 0.62$), higher-order topological lags (`nblag`), row-standardization ($W$), and inverse-distance gravity weights. |
| [**`Spatial-Autocorrelation.ipynb`**](file:///home/olande/Documents/GitHub/Spatial/Spatial-Autocorrelation.ipynb) | Global & Local Autocorrelation | Boston, MA (`boston_tracts.gpkg`, $n=506$) | Global Moran's $I$, Monte Carlo permutation distributions ($M=999$), Moran scatterplot regression equivalence, Local Moran ($I_i$), and Anselin LISA cluster maps (Hotspots vs. Coldspots). |
| [**`Bayesian-Spatial-Models.ipynb`**](file:///home/olande/Documents/GitHub/Spatial/Bayesian-Spatial-Models.ipynb) | Hierarchical Spatial Regression | Boston, MA (`boston_tracts.gpkg`, $n=506$) | Fast deterministic inference with INLA, classical BYM vs. modern scaled BYM2 models, spatial precision $\tau_s$, mixing parameter $\phi$, and non-linear dollar back-transformations (`inla.tmarginal`). |
| [**`Disease-Risk-Modeling.ipynb`**](file:///home/olande/Documents/GitHub/Spatial/Disease-Risk-Modeling.ipynb) | Small-Area Epidemiological Risk | Pennsylvania (`pennLC`, $n=67$ counties) | Demographic indirect standardization (16 strata), expected counts ($E_i$), raw Standardized Mortality Ratios (SMR), Poisson log-linear spatial models, Bayesian shrinkage, and exceedance probabilities ($P(RR > 1.20)$). |
| [**`Areal-Data-Issues.ipynb`**](file:///home/olande/Documents/GitHub/Spatial/Areal-Data-Issues.ipynb) | Methodological Hazards & Solutions | 256 simulated micro-units & regional strata | The Modifiable Areal Unit Problem (MAUP scale and zoning effects), the Ecological Fallacy (Robinson's paradox), spatial change of support (COSP), Bayesian melding, and downscaling. |

---

## Technical Standards

* **Package Management**: Managed strictly through `librarian::shelf(...)` (`spData`, `sf`, `spdep`, `INLA`, `SpatialEpi`, `tidyverse`, `patchwork`).
* **Visualization Philosophy**: Fully vectorized `ggplot2` multi-panel figures paired with `patchwork`. All cells separate variable assignments from visualization calls.
* **Canvas Geometry**: Tuned with explicit dual-kernel device dimensions (`ark.plot.*` for Positron and `repr.plot.*` for IRkernel/Jupyter) to guarantee snug, zero-margin rendering without scrollbars.
* **Documentation**: Independent, publication-grade mathematical and statistical explanations without external book dependencies.
