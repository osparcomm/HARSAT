# Checks convergence of an assessment model

Utilty function for use within assess_lmm. Checks whether a model has
converged. Currently only checks assessments with normal (or lognormal)
errors, where it considers whether:

- fixed effect estimates are away from their bounds

- random effect estimates are lower than their upper bounds; they can of
  course be equal to zero

- standard errors are present for model predictions and fixed effects
  estimates

- standard errors of the fixed effects estimates are not unrealistically
  small; the tolerance is chosen to be much smaller than seen in typical
  OSPAR assessments which have converged

  Model fits based on other distributions are assumed to have converged.
  Some checking is needed here in future

## Usage

``` r
check_convergence_lmm(assessment, coeff_se_tol = 0.001)
```

## Arguments

- assessment:

  An assessment from assess_lmm

- coeff_se_tol:

  The tolerance for checking whether standard errors on the fixed
  effects estimates are unrealistically small; defaults to 0.001

## Value

An integer: 0 indicates convergence, 1 indicates an issue
