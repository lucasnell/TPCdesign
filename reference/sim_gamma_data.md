# Simulate data with gamma distributed error

Simulate data with gamma distributed error

## Usage

``` r
sim_gamma_data(temp, n_reps, obs_cv, ctmin, ctmax, a, b, scale_tpc = FALSE)
```

## Arguments

- temp:

  Numeric vector giving the observed temperatures for each performance
  value.

- n_reps:

  Single integer for the number of experimental reps per temperature.

- obs_cv:

  Single numeric for the observation error coefficient of variation.
  Must be \>= 0. Observation error is simulated using a gamma
  distribution with mean equal to the expected performance. Where
  expected performance is zero (temperatures at or outside `ctmin` and
  `ctmax`, or `a = 0`), the simulated value is 0.

- ctmin:

  Single numeric for parameter `ctmin`.

- ctmax:

  Single numeric for parameter `ctmax`.

- a:

  Single numeric for parameter `a`.

- b:

  Single numeric for parameter `b`.

- scale_tpc:

  Single logical for whether to scale to make max value 1. Defaults to
  `FALSE`.

## Value

A dataframe of temperatures (column `"temp"`) and performance values
(column `"y"`) with gamma error.

## Examples

``` r
set.seed(1)
sim_gamma_data(temp = seq(10, 35, 5), n_reps = 3, obs_cv = 0.1,
               ctmin = 8, ctmax = 38, a = 1, b = 0.2)
#>    temp          y
#> 1    10   35.57286
#> 2    10   40.37871
#> 3    10   32.99959
#> 4    15  205.26965
#> 5    15  210.34050
#> 6    15  207.04828
#> 7    20  412.75265
#> 8    20  392.27028
#> 9    20  399.58658
#> 10   25  685.29039
#> 11   25  703.14054
#> 12   25  644.57431
#> 13   30 1079.00482
#> 14   30 1055.51676
#> 15   30 1089.18505
#> 16   35 1264.97387
#> 17   35 1180.10196
#> 18   35  949.37339
```
