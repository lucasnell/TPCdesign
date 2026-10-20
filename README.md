
<!-- README.md is generated from README.Rmd. Please edit that file -->

<!-- badges: start -->

[![R-CMD-check](https://github.com/lucasnell/TPCdesign/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/lucasnell/TPCdesign/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/lucasnell/TPCdesign/branch/main/graph/badge.svg)](https://app.codecov.io/gh/lucasnell/TPCdesign?branch=main)
<!-- badges: end -->

# TPCdesign

## Overview

Fitting a thermal performance curve (TPC) requires choosing which
temperatures to measure performance at. Temperatures placed in
uninformative parts of the curve can leave the critical thermal minimum
and maximum poorly estimated, even with many replicates.

`TPCdesign` chooses experimental temperatures that are expected to be
most informative for fitting a Brière-2 TPC, given rough prior guesses
at the curve’s parameters and how uncertain you are about them. The
optimization uses approximate coordinate exchange from the
[`acebayes`](https://cran.r-project.org/package=acebayes) package.

Main functions:

- `design_temps()`: choose informative temperatures.
- `briere2_tpc()` and `briere2_tpc_Topt()`: evaluate a Brière-2 curve
  and its thermal optimum.
- `sim_gamma_data()`: simulate performance data at a set of
  temperatures.

## Installation

`TPCdesign` contains C++ code, so installing from source requires a
compiler (Rtools on Windows, Xcode Command Line Tools on macOS).

``` r
# install.packages("remotes")
remotes::install_github("lucasnell/TPCdesign")
```

To also build the vignette:

``` r
remotes::install_github("lucasnell/TPCdesign", build_vignettes = TRUE)
```

### From source

Clone the repository:

``` bash
git clone https://github.com/lucasnell/TPCdesign.git
cd TPCdesign
```

Then, from R with the working directory set to the cloned folder,
install the dependencies and the package:

``` r
# install.packages("remotes")
remotes::install_deps(dependencies = TRUE)
remotes::install_local(".", build_vignettes = TRUE)
```

Or, from the terminal inside the cloned folder (requires the
dependencies listed in `DESCRIPTION` to already be installed):

``` bash
R CMD build .
R CMD INSTALL TPCdesign_*.tar.gz
```

## Usage

Give `design_temps()` your best guesses at `ctmin`, `ctmax`, `a`, and
`b`, plus how uncertain you are about `ctmin`, `ctmax`, and `log(b)`:

``` r
library(TPCdesign)

set.seed(1)
temps <- design_temps(
    n_temps = 8,
    ctmin = 8, ctmax = 38, a = 1, b = 0.2,
    ctmin_err = 2, ctmax_err = 1.5, logb_err = 0.26
)
temps
#> [1] 10.22 12.23 18.94 25.64 28.31 34.54 36.55 38.56
```

<img src="man/figures/README-design-plot-1.png" alt="" width="75%" />

Simulate data at those temperatures before running the experiment:

``` r
dat <- sim_gamma_data(temps, n_reps = 6L, obs_cv = 0.15,
                      ctmin = 8, ctmax = 38, a = 1, b = 0.2)
head(dat)
#>    temp        y
#> 1 10.22 43.34347
#> 2 10.22 38.28088
#> 3 10.22 41.01239
#> 4 10.22 41.54545
#> 5 10.22 46.95934
#> 6 10.22 48.94567
```

### Parallel processing

`design_temps()` runs its independent optimization starts (`n_starts`)
in parallel if you set a non-sequential
[`future`](https://cran.r-project.org/package=future) plan:

``` r
library(future)
plan(multisession, workers = 4L)
temps <- design_temps(n_temps = 8, ctmin = 8, ctmax = 38, a = 1, b = 0.2)
plan(sequential)
```

## Learn more

Full documentation is at <https://lucasnell.github.io/TPCdesign/>.

The vignette walks through the full workflow, including fitting
simulated data from an optimized and an evenly spaced design and
comparing parameter recovery:

``` r
vignette("TPCdesign", package = "TPCdesign")
```
