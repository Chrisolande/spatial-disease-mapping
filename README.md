# Spatial Data Science and Disease Mapping in R

R notebooks covering areal spatial data analysis, spatial autocorrelation, and Bayesian disease mapping.

---

## Notebooks

| Notebook | Focus | Data | What It Covers |
| :--- | :--- | :--- | :--- |
| [**`main.ipynb`**](file:///home/olande/Documents/GitHub/Spatial/main.ipynb) | US County Cancer Mortality ESDA | Contiguous US counties ($n=3,022$) from SEER and ACS | FIPS code matching, rate distribution, small-area variance (funnel plot), correlation and collinearity checks (VIF), Queen neighbor weights, Global Moran's $I$, OLS residual autocorrelation, Lagrange Multiplier error vs. lag tests, 6-order spatial correlogram, and LISA cluster maps. |
| [**`Spatial-Neighborhood-Matrices.ipynb`**](file:///home/olande/Documents/GitHub/Spatial/Spatial-Neighborhood-Matrices.ipynb) | Neighbor Weights & Graph Topology | Columbus, OH crime data ($n=49$) | Queen and Rook contiguity, $k$-nearest neighbors, distance bands, higher-order lags, and row-standardized weights ($W$). |
| [**`Spatial-Autocorrelation.ipynb`**](file:///home/olande/Documents/GitHub/Spatial/Spatial-Autocorrelation.ipynb) | Autocorrelation Diagnostics | Boston housing tracts ($n=506$) | Global Moran's $I$, Monte Carlo permutation tests ($M=999$), Moran scatterplot, and Local Moran (LISA) hotspot and coldspot maps. |
| [**`Bayesian-Spatial-Models.ipynb`**](file:///home/olande/Documents/GitHub/Spatial/Bayesian-Spatial-Models.ipynb) | Spatial Regression with INLA | Boston housing tracts ($n=506$) | Fitting Besag-York-Mollié (BYM and BYM2) models with R-INLA, spatial precision parameters, and posterior marginals. |
| [**`Disease-Risk-Modeling.ipynb`**](file:///home/olande/Documents/GitHub/Spatial/Disease-Risk-Modeling.ipynb) | Small-Area Disease Risk | Pennsylvania lung cancer data ($n=67$) | Indirect age/race/gender rate standardization, expected counts ($E_i$), raw SMRs, Poisson spatial regression, and exceedance probabilities ($P(RR > 1.20)$). |
| [**`Areal-Data-Issues.ipynb`**](file:///home/olande/Documents/GitHub/Spatial/Areal-Data-Issues.ipynb) | Common Spatial Pitfalls | Simulated grid data | Modifiable Areal Unit Problem (MAUP scale and zoning effects), ecological fallacy, and spatial misalignment. |

---

## Data

The `data/` folder contains:

* **`cancer_reg.csv`**: County-level age-adjusted cancer mortality and socioeconomic predictors from SEER and the American Community Survey.
* **`death.csv` & `incd.csv`**: Reference tables used to match county names to 5-digit FIPS codes.
* **`us_counties.graph`**: Adjacency graph file for the 3,022 contiguous US counties, formatted for `inla.read.graph()`.

---

## Requirements

Packages can be loaded using `librarian` (which installs anything missing):

```r
if (!requireNamespace("librarian", quietly = TRUE)) install.packages("librarian")
librarian::shelf(sf, spdep, tigris, tidyverse, patchwork, viridis, SpatialEpi)
```

To run the Bayesian models in `Bayesian-Spatial-Models.ipynb` and `Disease-Risk-Modeling.ipynb`, install R-INLA:

```r
install.packages("INLA", repos = c(getOption("repos"), INLA = "https://inla.r-inla-download.org/R/stable"))
```
