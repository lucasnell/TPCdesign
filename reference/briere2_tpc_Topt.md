# Thermal optimum for Brière-2 thermal performance curve (TPC)

Note that parameter `a` is not necessary for this calculation.

## Usage

``` r
briere2_tpc_Topt(ctmin, ctmax, b)
```

## Arguments

- ctmin:

  Numeric vector for parameter `ctmin`.

- ctmax:

  Numeric vector for parameter `ctmax`.

- b:

  Numeric vector for parameter `b`.

## Value

A numeric vector of optimum temperatures.

## Examples

``` r
briere2_tpc_Topt(ctmin = 8, ctmax = 38, b = 0.2)
#> [1] 34.95608
```
