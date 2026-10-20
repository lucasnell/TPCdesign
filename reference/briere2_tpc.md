# Brière-2 thermal performance curve (TPC)

Note that this TPC does not perform well when ctmin is positive, but
there are negative temperatures. For example, if you run
`briere2_tpc(c(-10, 0, 10), ctmin = 5, ctmax = 30, a = 1, b = 0.2, TRUE)`,
you'll see that a temperature of `-10` gives a greater performance value
than `10`.

## Usage

``` r
briere2_tpc(temp, ctmin, ctmax, a, b, scale = FALSE)
```

## Arguments

- temp:

  Numeric vector of temperatures

- ctmin:

  Single numeric for parameter `ctmin`.

- ctmax:

  Single numeric for parameter `ctmax`.

- a:

  Single numeric for parameter `a`.

- b:

  Single numeric for parameter `b`.

- scale:

  Single logical for whether to scale to make max value 1. Defaults to
  `FALSE`.

## Value

A numeric vector for measure of performance for each in `temp`

## Examples

``` r
briere2_tpc(temp = seq(10, 35, 5), ctmin = 8, ctmax = 38, a = 1, b = 0.2)
#> [1]   38.94589  196.57798  427.82459  709.86800 1000.37293 1177.21574
```
