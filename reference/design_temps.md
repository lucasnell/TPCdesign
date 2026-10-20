# Design experimental temperatures

Design experimental temperatures

## Usage

``` r
design_temps(
  n_temps,
  ctmin,
  ctmax,
  a,
  b,
  ctmin_err = 2.5,
  ctmax_err = 1.5,
  logb_err = 0.26,
  min_sep = 2,
  n_filler = 1L,
  n_draws = 250L,
  n_starts = 7L,
  digits = 2L,
  min_temp = NULL,
  max_temp = NULL
)
```

## Arguments

- n_temps:

  Single integer for the total number of temperatures to return,
  including any filler temperatures added via `n_filler`. Must be \>= 2.

- ctmin:

  Single numeric for your best guess at the curve's critical thermal
  minimum, i.e., the lower temperature at which performance reaches
  zero. Used as the center of the uniform prior for `ctmin` (see
  `ctmin_err`).

- ctmax:

  Single numeric for your best guess at the curve's critical thermal
  maximum, i.e., the upper temperature at which performance reaches
  zero. Used as the center of the uniform prior for `ctmax` (see
  `ctmax_err`). Must be \> `ctmin`.

- a:

  Single numeric for your best guess at the Brière-2 curve's `a`
  parameter, which scales the overall magnitude of performance. Unlike
  `ctmin`, `ctmax`, and `b`, `a` is treated as known exactly (there is
  no corresponding `*_err` argument) since the shape of the optimal
  design is insensitive to it. Must be \> 0.

- b:

  Single numeric for your best guess at the Brière-2 curve's `b`
  parameter, which controls the asymmetry (skew) of the curve between
  `ctmin` and `ctmax`. Used as the center of the uniform prior for
  `log(b)` (see `logb_err`). Must be \> 0.

- ctmin_err:

  Single numeric for how uncertain you are about `ctmin`. This means
  priors for `ctmin` will be generated from a uniform distribution with
  minimum `ctmin - ctmin_err` and maximum `ctmin + ctmin_err`. Larger
  values reflect less confidence in your `ctmin` guess and cause the
  design to hedge across a wider range of possible curves. Must be
  `>= 0`. Defaults to `2.50`.

- ctmax_err:

  Single numeric for how uncertain you are about `ctmax`. This means
  priors for `ctmax` will be generated from a uniform distribution with
  minimum `ctmax - ctmax_err` and maximum `ctmax + ctmax_err`. Larger
  values reflect less confidence in your `ctmax` guess and cause the
  design to hedge across a wider range of possible curves. Must be
  `>= 0`. Defaults to `1.50`.

- logb_err:

  Single numeric for how uncertain you are about `b`, expressed on the
  log scale. This means priors for `log(b)` will be generated from a
  uniform distribution with minimum `log(b) - logb_err` and maximum
  `log(b) + logb_err`. Larger values reflect less confidence in your `b`
  guess and cause the design to hedge across a wider range of possible
  curve shapes. Must be `>= 0`. Defaults to `0.26`.

- min_sep:

  Single numeric specifying the minimum separation (in temperature
  units) enforced between optimized temperatures, to avoid wasting
  experimental effort on temperatures that are too close together to be
  useful. The step where equally-spaced temperatures are added (only
  happens when `n_filler >= 1`) is not affected by this argument, and
  unless the number of points is very high relative to the difference in
  maximum (`ctmax + ctmax_err`) and minimum (`ctmin + ctmin_err`)
  temperatures surveyed, having a minimum separation shouldn't change
  this step anyway. Defaults to `2`.

- n_filler:

  Single integer specifying how many of the `n_temps` total temperatures
  should be equally spaced (chosen to fill the largest gaps left by the
  optimized temperatures) rather than individually optimized. Equally
  spacing points can be a good bet hedging strategy if you're quite
  unsure of the TPC parameters (`ctmin`, `ctmax`, `a`, `b`) you're
  using, since it guards against the optimization overfitting to a
  misspecified curve. Must be between `0` and `n_temps`. Defaults to
  `1L`.

- n_draws:

  Single integer specifying the number of Monte Carlo draws taken from
  the priors on `ctmin`, `ctmax`, and `b` to approximate the expected
  information gained by a candidate design. Increase if your certainty
  is low (i.e., `*_err` parameters are high), since wider priors need
  more draws to be well approximated. Larger values give a more accurate
  (but slower) optimization. Defaults to `250L`.

- n_starts:

  Single integer specifying the number of independent random starting
  designs the coordinate-exchange optimization is run from, with the
  best-performing result across all starts returned. Increasing this
  reduces the chance the optimizer settles on a locally (rather than
  globally) optimal set of temperatures, at the cost of more
  computation; this cost can be mitigated by running starts in parallel
  (see Examples). Defaults to `7L`.

- digits:

  Single integer specifying the number of decimal places to round the
  returned temperatures to. Defaults to `2L`.

- min_temp:

  Single numeric specifying the minimum possible value of design
  temperatures to allow. If `NULL` (the default), this is set to
  `ctmin - ctmin_err`. Must be `<= ctmin`.

- max_temp:

  Single numeric specifying the maximum possible value of design
  temperatures to allow. If `NULL` (the default), this is set to
  `ctmax + ctmax_err`. Must be `>= ctmax`.

## Value

A vector of temperatures at which to sample.

## Examples

``` r
design_temps(n_temps = 10, ctmin = 5, ctmax = 40, a = 1, b = 0.2)
#>  [1] 10.30 12.31 14.32 20.80 27.28 29.29 34.46 36.47 38.48 40.56

# To run in parallel:
# library(future.apply)
# plan(multisession, workers = 4L)
# design_temps(n_temps = 10, ctmin = 5, ctmax = 40, a = 1, b = 0.2)
```
