# ctsm.VDS.loglik.calc

ctsm.VDS.loglik.calc

## Usage

``` r
ctsm.VDS.loglik.calc(
  theta,
  data,
  index.theta,
  minus.twice = FALSE,
  cumulate = FALSE
)
```

## Arguments

- theta:

  The maximum imposex stage

- data:

  Individual imposex data.

- index.theta:

  Allows for stage names that do not run from 0 to `theta` (optional)

- minus.twice:

  A logical specifying whehter to calculated the liklihood (FALSE) or
  the deviance (TRUE); default FALSE

- cumulate:

  A logical specifying whether to use cumulative probabilities; default
  FALSE
