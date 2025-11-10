# Default control parameters for `run_assessment`

**\[experimental\]** Default parameters that control the way the
assessment is run. Presently only includes parameters for post-hoc
power, but it is intended to move `recent_trend` here, along with the
arguments that control the calculation of numerical derivatives.

## Usage

``` r
run_control_default()
```

## Value

A list with the following components:

- `power` A list with the following components (all expressed as
  percentages):

  - `target_power` default = 90%

  - `target_trend` default = 5%

  - `size` default = 5% The power calculations are currently only
    applied to log-normally distributed data, which is why the trend is
    expressed as a percentage.

- `recent_change` A list with the following components:

  - `n_year_fit` default = 5L

  - `n_year_positive` default = 5L A recent change will only be computed
    if there are at least `n_year_fit` years of data in the recent
    period, of which at least `n_year_positive` contain at least one
    non-censored measurement. This only affects normally or log-normally
    distributed data.
