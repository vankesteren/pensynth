
# Penalized synthetic control estimation

[![R-CMD-check](https://github.com/vankesteren/pensynth/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/vankesteren/pensynth/actions/workflows/R-CMD-check.yaml)

The goal of `pensynth` is to make it easier to perform penalized synthetic control in the spirit of [Abadie & L'Hour (2021)](https://doi.org/10.1080/01621459.2021.1971535). 

## Features
- Faster than the original `Synth::synth` implementation for "vanilla" synthetic controls, even for large donor pools.
- Built-in cross-validation on the pre-intervention outcome timeseries to determine the penalty parameter. (see example).
- Plotting of the full solution path for cross-validated penalized synthetic controls.

NB: in this implementation, variable weights have to be pre-specified (unlike in the original synthetic control implementaion). Additionally, currently only a single treated unit is supported. 

## Installation

You can install the development version of `pensynth` like so:

``` r
remotes::install_github("vankesteren/pensynth")
```

## Example

``` r
library(pensynth)
set.seed(45)

# Generate some data
N_covar <- 7   # number of covariates
N_donor <- 50  # number of donor units
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

# References

Abadie, A., & L’Hour, J. (2021). A penalized synthetic control estimator for disaggregated data. _Journal of the American Statistical Association, 116_(536), 1817-1834.

Some of the code was inspired by (but heavily adapted from) [jeremylhour/pensynth](https://github.com/jeremylhour/pensynth). Where this was the case, the code is commented.
