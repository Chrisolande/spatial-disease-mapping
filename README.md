# Spatial Data Science & Bayesian Disease Mapping in R

[![R Version](https://img.shields.io/badge/R-%E2%89%A54.6-276DC3?logo=r&logoColor=white)](https://www.r-project.org/)
[![Spatial Stack](https://img.shields.io/badge/Stack-sf%20%7C%20spdep%20%7C%20tigris-059669)](https://r-spatial.github.io/sf/)
[![Bayesian Engine](https://img.shields.io/badge/Bayesian-R--INLA-e11d48)](https://www.r-inla.org/)
[![License](https://img.shields.io/badge/License-MIT-slate)](LICENSE)

A collection of notebooks covering areal spatial statistics, neighbor network topology, spatial autocorrelation diagnostics, and hierarchical Bayesian disease mapping in R. 

The repository includes both modular methodological reference notebooks and a full-scale empirical study (`main.ipynb`) evaluating age-adjusted cancer mortality across $n = 3,022$ contiguous US counties using clinical registries from NCI/CDC SEER and demographic covariates from the US Census American Community Survey (ACS).

---

## Repository Structure

```text
.
├── main.ipynb                           # Nationwide US county cancer mortality ESDA
├── Spatial-Neighborhood-Matrices.ipynb  # Contiguity weights, distance bands, and graph metrics
├── Spatial-Autocorrelation.ipynb        # Global Moran's I, permutation tests, and LISA clusters
├── Bayesian-Spatial-Models.ipynb        # Hierarchical spatial regression with R-INLA (BYM/BYM2)
├── Disease-Risk-Modeling.ipynb          # Epidemiological rate standardization and Poisson disease mapping
├── Areal-Data-Issues.ipynb             # Modifiable Areal Unit Problem (MAUP) and ecological fallacy
└── data/
    ├── cancer_reg.csv                   # County cancer mortality rates and 33 ACS covariates
    ├── death.csv                        # SEER mortality reference table with 5-digit FIPS codes
    ├── incd.csv                         # SEER incidence reference table with 5-digit FIPS codes
    ├── cancer_data_notes.csv            # Variable definitions and registry documentation
    └── us_counties.graph                # Binary topological adjacency graph for 3,022 counties
```

---

## Notebook Guide

| Notebook | Focus | Dataset | Description |
| :--- | :--- | :--- | :--- |
| [**`main.ipynb`**](main.ipynb) | Nationwide Cancer Mortality ESDA | Contiguous US counties ($n=3,022$) | Full exploratory pipeline: FIPS matching, age-standardization check, distribution moments, small-area variance instability (funnel plot), count overdispersion ($\text{Var}/\text{E} \approx 1,553$), bivariate screening under Bonferroni/Holm FWER control, collinearity diagnostics (VIF), Queen contiguity graph generation, Global Moran's $I$ ($0.462$), OLS residual autocorrelation ($0.342, p < 10^{-200}$), Lagrange Multiplier error vs. lag tests, 6-order spatial correlogram, and Local Moran (LISA) cluster mapping with FDR and Bonferroni error control. |
| [**`Spatial-Neighborhood-Matrices.ipynb`**](Spatial-Neighborhood-Matrices.ipynb) | Spatial Adjacency & Topology | Columbus, OH crime data ($n=49$) | Constructing and comparing Queen and Rook polygon contiguity, $k$-nearest neighbors, distance bands, higher-order topological lags (`nblag`), and row-standardized weight matrices ($W$). |
| [**`Spatial-Autocorrelation.ipynb`**](Spatial-Autocorrelation.ipynb) | Autocorrelation Diagnostics | Boston housing tracts ($n=506$) | Global Moran's $I$, Monte Carlo permutation distributions ($M=999$), regression slope equivalence on Moran scatterplots, and Anselin Local Moran ($I_i$) hotspot/coldspot classification. |
| [**`Bayesian-Spatial-Models.ipynb`**](Bayesian-Spatial-Models.ipynb) | Hierarchical Spatial Regression | Boston housing tracts ($n=506$) | Fitting Besag-York-Mollié (BYM and scaled BYM2) models with R-INLA, estimating spatial precision ($\tau_s$) and mixing parameter ($\phi$), and extracting posterior marginal distributions. |
| [**`Disease-Risk-Modeling.ipynb`**](Disease-Risk-Modeling.ipynb) | Small-Area Disease Risk | Pennsylvania lung cancer data ($n=67$) | Indirect age, gender, and race rate standardization across 16 demographic strata, computing expected counts ($E_i$), raw SMRs, Poisson log-linear models, Bayesian shrinkage, and exceedance probabilities ($P(RR > 1.20)$). |
| [**`Areal-Data-Issues.ipynb`**](Areal-Data-Issues.ipynb) | Spatial Aggregation Hazards | 256 simulated micro-units | Demonstrating the Modifiable Areal Unit Problem (MAUP scale and zoning effects), ecological fallacy (Robinson's paradox), spatial change of support, and downscaling principles. |

---

## Mathematical Formulations

### 1. Global Moran's $I$ (Spatial Autocorrelation)
Measures the overall linear association between an observed continuous variable $y$ and its spatially lagged neighbor values $W y$:

$$I = \frac{n}{\sum_{i=1}^n \sum_{j=1}^n w_{ij}} \frac{\sum_{i=1}^n \sum_{j=1}^n w_{ij}(y_i - \bar{y})(y_j - \bar{y})}{\sum_{i=1}^n (y_i - \bar{y})^2}$$

where $w_{ij}$ are elements of a row-standardized spatial weights matrix $W$ ($w_{ij} > 0$ if counties $i$ and $j$ share a boundary, and $\sum_j w_{ij} = 1$). Under the null hypothesis of spatial randomness ($H_0$), the expected value is $E[I] = -1 / (n - 1) \approx -0.0003$.

### 2. Anselin Local Moran ($I_i$)
Decomposes global spatial autocorrelation into local contributions to detect localized spatial clusters and spatial outliers:

$$I_i = \frac{z_i - \bar{z}}{s^2} \sum_{j=1}^n w_{ij} (z_j - \bar{z})$$

Counties are categorized into four quadrants:
* **High-High (Hotspot):** Elevated rate surrounded by elevated neighbor rates.
* **Low-Low (Coldspot):** Low rate surrounded by low neighbor rates.
* **High-Low / Low-High (Spatial Outliers):** A county whose rate diverges significantly from its neighbors.

### 3. Hierarchical BYM2 Model (Riebler et al., 2016)
For disease mapping across areal regions, the scaled Besag-York-Mollié (BYM2) model decomposes area-level risk into structured spatial clustering and unstructured noise:

$$y_i = \mathbf{x}_i^T \boldsymbol{\beta} + \frac{1}{\sqrt{\tau}} \left(\sqrt{1 - \phi} \, v_i + \sqrt{\phi} \, u_i^*\right)$$

where:
* $\tau$: Overall marginal precision parameter.
* $\phi \in [0, 1]$: Spatial mixing parameter measuring the proportion of marginal variance explained by the spatial network ($u_i^*$) versus unstructured noise ($v_i$).
* $u_i^*$: Scaled Intrinsic Autoregressive (ICAR) component:
  $$u_i \mid u_{-i} \sim \mathcal{N}\left(\frac{\sum_{j \in \mathcal{N}(i)} u_j}{d_i}, \frac{1}{\tau_u d_i}\right)$$
  scaled such that the geometric mean of its marginal variances equals 1 ($\frac{1}{n}\sum_{i=1}^n \text{Var}(u_i^*) \approx 1$).
* $v_i \sim \mathcal{N}(0, 1)$: Independent and identically distributed Gaussian white noise.

### 4. Penalized Complexity (PC) Priors (Simpson et al., 2017)
Priors on the precision and mixing parameters are assigned using Penalized Complexity (PC) formulations that decay exponentially with distance from a base model:
* **Precision ($\tau$):** $\text{Prob}(1/\sqrt{\tau} > 0.5) = 0.01$ (penalizes excessive random effect variance).
* **Mixing parameter ($\phi$):** $\text{Prob}(\phi < 0.5) = 0.5$ (neutral, uninformative prior regarding structured vs. unstructured dominance).

---

## Key Empirical Findings from `main.ipynb`

* **Distribution Properties:** Age-adjusted cancer mortality centers symmetrically at a mean of $178.66$ and median of $178.10$ per 100k, with a modest skewness of $0.27$ and an IQR of $34.0$. 
* **Small-Area Noise:** Sparser counties ($< 20,000$ residents) show high crude rate dispersion ($80$ to $360$ per 100k) due to small population denominators, demonstrating the need for Bayesian shrinkage.
* **Predictor Collinearity:** Median household income and poverty rate correlate at $r = -0.79$, resulting in Variance Inflation Factors of $\text{VIF} \approx 5.05$ (poverty) and $\text{VIF} \approx 4.74$ (income), with a design matrix condition number of $\kappa = 6.30$.
* **Persistent Spatial Autocorrelation:** Raw county mortality exhibits strong spatial clustering ($I = 0.462, z = 41.95, p < 2.2 \times 10^{-16}$). In an OLS regression conditioned on poverty, income, education, and health coverage, the residuals retain significant autocorrelation ($\text{Residual Moran's } I = 0.342, z = 31.25, p < 10^{-200}$).
* **Spatial Error Process:** Lagrange Multiplier diagnostics confirm that unobserved spatial confounders dominate contagious diffusion ($\text{Robust LM Error} = 104.51 \gg \text{Robust LM Lag} = 42.32, p < 2.2 \times 10^{-16}$), justifying an areal spatial error / random effects model.
* **Spatial Range:** Autocorrelation decays gradually across contiguity orders ($I_1 = 0.462 \to I_6 = 0.286$, spanning roughly $400\text{ km}$), indicating regional clustering across Appalachia and the Mississippi Delta.

---

## Setup & Reproduction

### Prerequisites

All notebooks run in R 4.6+ using either Positron (Ark kernel) or Jupyter.

### Package Installation

Libraries can be loaded and auto-installed via `librarian`:

```r
if (!requireNamespace("librarian", quietly = TRUE)) {
  install.packages("librarian")
}

librarian::shelf(
  sf, spdep, tigris, tidyverse, patchwork, viridis, SpatialEpi,
  quiet = TRUE
)
```

To run the Bayesian models in `Bayesian-Spatial-Models.ipynb` and `Disease-Risk-Modeling.ipynb`, install R-INLA:

```r
install.packages(
  "INLA",
  repos = c(getOption("repos"), INLA = "https://inla.r-inla-download.org/R/stable")
)
```

### Running the Notebooks

1. Clone the repository:
   ```bash
   git clone https://github.com/Chrisolande/spatial-disease-mapping.git
   cd spatial-disease-mapping
   ```
2. Open in Positron, RStudio, or launch via Jupyter:
   ```bash
   jupyter lab main.ipynb
   ```

---

## References

* Besag, J., York, J., & Mollié, A. (1991). Bayesian image restoration, with two applications in spatial statistics. *Annals of the Institute of Statistical Mathematics*, 43(1), 1–20.
* Riebler, A., Sørbye, S. H., Simpson, D., & Rue, H. (2016). An intuitive Bayesian spatial model for disease mapping that accounts for scaling. *Statistical Methods in Medical Research*, 25(4), 1145–1165.
* Simpson, D., Rue, H., Riebler, A., Martins, T. G., & Sørbye, S. H. (2017). Penalising model component complexity: A principled, practical approach to constructing priors. *Statistical Science*, 32(1), 1–28.
* Wu, K. L., & Banerjee, S. (2025). Spatial Confounding in Multivariate Areal Data Analysis. *Biometrics*, arXiv:2505.07232.
* Anselin, L. (1995). Local Indicators of Spatial Association—LISA. *Geographical Analysis*, 27(2), 93–115.
