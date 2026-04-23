# aftPenCDA

`aftPenCDA` is an R package for fitting penalized accelerated failure time (AFT) models using induced smoothing and coordinate descent algorithms. Computationally intensive components are implemented in 'C++' via 'Rcpp' (RcppArmadillo backend) to ensure scalability in high-dimensional settings.

The package supports both right-censored survival data and clustered partly interval-censored survival data, and provides flexible variable selection through several penalty functions.

---

## Features

- Induced smoothing-based AFT estimation for nonsmooth estimating equations
- Penalized coordinate descent algorithm
- Efficient C++ implementation via Rcpp (RcppArmadillo)
- Support for high-dimensional regression
- Supported penalties:
  - Broken Adaptive Ridge (BAR)
  - LASSO
  - Adaptive LASSO (ALASSO)
  - SCAD
- Variance estimation methods:
  - `"CF"`: closed-form plug-in estimator
  - `"ZL"`: perturbation-based estimator based on Zeng and Lin (2008)

---

## Installation

You can install the development version from GitHub:

```r
# install.packages("devtools")
devtools::install_github("seonsy/aftPenCDA")
```

## Main functions

### `aftpen()`

Fits a penalized AFT model for right-censored survival data.

### `aftpen_pic()`

Fits a penalized AFT model for clustered partly interval-censored survival data.

Both functions share the same interface:

```r
aftpen(dt, lambda = 0.1, se = "CF", type = "BAR")
aftpen_pic(dt, lambda = 0.1, se = "CF", type = "BAR")
```
## Input data format

### Right-censored data (`aftpen()`)

- `dt`: a data frame where:
  - first column: `y` (observed time)
  - second column: `d` (event indicator; `1 = event`, `0 = censoring`)
  - remaining columns: covariates

---

### Clustered partly interval-censored data (`aftpen_pic()`)

- `dt`: a data frame containing:
  - `L`, `R`: interval endpoints
  - `delta`: exact observation indicator (`1 = exact`, `0 = censored`)
  - `id`: cluster identifier
  - remaining columns: covariates

## Algorithm

The method combines induced smoothing with a coordinate descent algorithm. 
A quadratic approximation is constructed via Cholesky decomposition, leading to a least-squares-type problem.

Efficient coordinate-wise updates are then applied under different penalties.

## Example

### Right-censored data
```r
library(aftPenCDA)

data("simdat_rc")

fit <- aftpen(simdat_rc, lambda = 0.1, se = "CF", type = "BAR")

fit$beta
```

### Clustered partly interval-censored data
```r
data("simdat_pic")

fit_pic <- aftpen_pic(simdat_pic, lambda = 0.001, se = "CF", type = "BAR")

fit_pic$beta
```

## Arguments

| Argument   | Description |
|-----------|------------|
| `lambda`   | Tuning parameter controlling penalization strength |
| `type`     | `"BAR"`, `"LASSO"`, `"ALASSO"`, `"SCAD"` |
| `se`       | Variance estimation method (`"CF"` or `"ZL"`) |
| `r`        | SCAD tuning parameter (default: `3.7`) |
| `eps`      | Convergence tolerance (default: `1e-8`) |
| `max.iter` | Maximum number of iterations (default: `100`) |

## Value

Both functions return a list with components:

- `beta`: final penalized coefficient estimate

## References

Wang, You-Gan, and Yudong Zhao. 2008. “Weighted Rank Regression for Clustered Data Analysis.” *Biometrics* 64 (1): 39–45.

Dai, L., K. Chen, Z. Sun, Z. Liu, and G. Li. 2018. “Broken Adaptive Ridge Regression and Its Asymptotic Properties.” *Journal of Multivariate Analysis* 168: 334–351.

Zeng, Donglin, and D. Y. Lin. 2008. “Efficient Resampling Methods for Nonsmooth Estimating Functions.” *Biostatistics* 9 (2): 355–363.



## Note

This package is under development. Functionality and interfaces may change in future versions.
