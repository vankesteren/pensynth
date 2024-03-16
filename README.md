
# Penalized synthetic control estimation

[![R-CMD-check](https://github.com/vankesteren/pensynth/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/vankesteren/pensynth/actions/workflows/R-CMD-check.yaml)
[![pensynth status badge](https://vankesteren.r-universe.dev/badges/pensynth)](https://vankesteren.r-universe.dev/pensynth)
[![cran version](
https://www.r-pkg.org/badges/version/pensynth)](https://cran.r-project.org/package=pensynth)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)

The goal of `pensynth` is to make it easier to perform penalized synthetic control in the spirit of [Abadie & L'Hour (2021)](https://doi.org/10.1080/01621459.2021.1971535). 

## Features
- Faster than the original `Synth::synth` implementation for "vanilla" synthetic controls, even for large donor pools, because we use the [`clarabel`](https://oxfordcontrol.github.io/ClarabelDocs/stable/) quadratic program solver.
- Built-in cross-validation on the pre-intervention outcome timeseries to determine the penalty parameter (see example below).
- Plotting of the full solution path for cross-validated penalized synthetic controls.

NB: in this implementation, variable weights have to be pre-specified (unlike in the original synthetic control implementation). Additionally, currently only a single treated unit is supported. 

## Installation

I recommend installing `pensynth` from r-universe like so:

```r
install.packages("pensynth", repos = c("https://vankesteren.r-universe.dev", "https://cloud.r-project.org"))
```

The package is also on CRAN, so it can be installed directly from there too.

```r
install.packages("pensynth")
```

The latest development version can also be installed directly from this repository:

```r
remotes::install_github("vankesteren/pensynth")
```

## Why penalization?

Penalized synthetic control yields a smooth transition between the synthetic control method (when $\lambda = 0$) and nearest neighbour matching (when $\lambda \to \infty$).

When the treated unit is in the convex hull of the donor units (which is more likely when there are many donors) there is no unique solution for the unit weights of synthetic control. In these cases especially, the penalty can help because it prefers solutions with donors closer in covariate space. 

The `pensynth` implementation achieves this through optimizing the following objective:

```math
\min_{\boldsymbol{w}} \left[ \| \boldsymbol{x}_1 - \boldsymbol{X}_0 \boldsymbol{w} \|^2 + \lambda \sum_{d\in D} \boldsymbol{w}_d \|\boldsymbol{x}_1 - \boldsymbol{x}_{d}\|^2 \right]
```

```math
\text{s.t.} \quad \boldsymbol{w}_1 \geq 0, ..., \boldsymbol{w}_D \geq 0,
\, \sum_{d\in D} \boldsymbol{w}_d = 1
```

Where 
- $\boldsymbol{x}_1$ is the column vector of treated unit covariates,  
- $\boldsymbol{X}_0$ are the covariate values for the donor units, 
- $D$ is the number of donor units, 
- $\lambda$ is the penalty parameter, and 
- $\boldsymbol{x}_{d}$ is taken to be the $d^{th}$ column of $\boldsymbol{X}_0$.

The first term in the objective is the same (up to variable weights) as the original synthetic control, and the second term is the nearest neighbour matching penalty.

## Example

``` r
library(pensynth)
set.seed(45)

# Generate some data
N_covar  <- 7   # number of covariates
N_donor  <- 50  # number of donor units
N_target <- 12 # pre-intervention time series length
w  <- runif(N_donor) # true unit weights
w[5:N_donor] <- 0    # set most to 0 (sparse)
w  <- w / sum(w)     # normalize unit weights
v  <- rep(1, N_covar) # equal variable weights

# Create the synthetic control data matrices
# following the Synth::synth() notation
X0 <- matrix(rnorm(N_covar*N_donor), N_covar)  
X1 <- X0 %*% w
Z0 <- matrix(rnorm(N_target*N_donor), N_target)
Z1 <- Z0 %*% w

# Run penalized synthetic control
# estimate lambda using pre-intervention timeseries MSE
res <- cv_pensynth(X1, X0, v, Z1, Z0)
plot_path(res)
```
![cvplot](img/cvplot.png)

```r
# Compare final weights to true weights
plot(w, main = "Estimated and true unit weights", xlab = "Donor unit", ylab = "Weight")
points(res$w_opt, pch = 3)
```
![wplot](img/weights.png)

# References

Abadie, A., & L’Hour, J. (2021). A penalized synthetic control estimator for disaggregated data. _Journal of the American Statistical Association, 116_(536), 1817-1834.

Some of the code was inspired by (but heavily adapted from) [jeremylhour/pensynth](https://github.com/jeremylhour/pensynth). Where this was the case, the code is commented.
