# Estimates cut-points in ordinal model (imposex assessments)

Estimates cut-points in ordinal model (imposex assessments)

## Usage

``` r
ordinal_theta_est(data, theta = NULL, ref_level = NULL, calc_vcov = FALSE)
```

## Arguments

- data:

  individual ordinal data (e.g. for imposex assessments)

- theta:

  an (optional) set of intial parameter values for theta; defaults to
  NULL

- ref_level:

  an (optional) reference level for parameter estimation; defaults to a
  index with intermediate levels of imposex

- calc_vcov:

  logical specifying whether to calculate the covariance matrix of the
  parameter estimates; defaults to FALSE
