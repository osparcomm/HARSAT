# Calculates confidence limits for imposex time series

Calculates confidence limits for imposex time series

## Usage

``` r
ordinal_theta_cl(fit, nsim = 1000, annual_indices = NULL)
```

## Arguments

- fit:

  The output from a call to `ordinal_theta_est` (sort of)

- nsim:

  The number of simulations on which each set of confidence limits is
  based; default 1000
