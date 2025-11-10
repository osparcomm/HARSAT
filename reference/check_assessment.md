# Check whether the assessments have converged

Checks whether the assessments in a HARSAT assessment object have
converged. Currently only does detailed checks for models with normal or
lognormal errors (all chemical timeseries and some biological effects).
Here it checks whether:

- fixed effect estimates are away from their bounds

- random effect estimates are lower than their upper bounds; they can of
  course be equal to zero

- standard errors are present for model predictions and fixed effects
  estimates

- standard errors of the fixed effects estimates are realistic (very
  small standard errors indicate problems with the numerical
  differencing used to compute the standard errors)

## Usage

``` r
check_assessment(assessment_ob, save_result = FALSE)
```

## Arguments

- assessment_ob:

  A HARSAT assessment object resulting from a call to ctsm_assessment

- save_result:

  Saves the identifiers of the timeseries that have not converged;
  defaults to FALSE. When an assessment has been done in stages, the
  output also identifies those timeseries that have not yet been
  assessed

## Value

A list of two character vectors:

- not_converged identifies timeseries that have not converged

- not_assessed identifies timeseries that have not been assessed
