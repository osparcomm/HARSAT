# Assess timeseries for trends and status

Fits a model to each timeseries, test for any temporal trend and compare
with thresholds. Need to add a lot more in details.

## Usage

``` r
run_assessment(
  ctsm_ob,
  subset = NULL,
  AC = NULL,
  get_AC_fn = NULL,
  recent_trend = 20L,
  parallel = FALSE,
  extra_data = NULL,
  control = list(),
  ...
)
```

## Arguments

- ctsm_ob:

  A HARSAT object resulting from a call to create_timeSeries

- subset:

  An optional vector specifying which timeseries are to be assessed.
  Might be used if the assessment is to be done in chunks because of
  size, or when refitting a timeseries model which has not converged. An
  expression will be evaluated in the timeSeries component of ctsm_ob;
  use 'series' to identify individual timeseries.

- AC:

  A character vector identifying the thresholds to be used in status
  assessments. These should be in the threshold reference table.
  Defaults to NULL; i.e. no thresholds are used.

- get_AC_fn:

  An optional function that overrides get_AC_default. See details (which
  need to be written).

- recent_trend:

  An integer giving the number of years which are used in the assessment
  of recent trends. For example, a value of 20 (the default) consider
  trends in the last twenty year.

- parallel:

  A logical which determines whether to use parallel computation;
  default = FALSE.

- extra_data:

  **\[experimental\]** A named list used to pass additional data to
  specific assessment routines. At present it is only used for imposex
  assessments, where it passes two data frames called `VDS_estimates`
  and `VDS_confidence_limits`. Defaults to NULL, This argument will be
  generalised in the near future, so expect it to change.

- control:

  **\[experimental\]** A list of control parameters that allow the user
  to modify the way the assessment is run. These currently include
  parameters for post-hoc power calculations and for the calculation of
  the 'recent_change' (in addition to the `recent_trend` argument
  above). See details.

- ...:

  Extra arguments which are passed to assessment_engine. See details
  (which need to be written).

## Details

### Control parameters

Some aspects of the model output are controlled using parameters which
can be modified using the `control` argument. The default parameter list
is described below. The `control` argument only needs to specify any
changes. For example, to change the target power from 90% (default) to
80%, use

`control = list(power = list(target_power = 80))`

The default parameters are a list with the following components:

`power`

A list with the following components (all expressed as percentages):

- `target_power` default = 90

- `target_trend` default = 5

- `size` default = 5

These affect the post-how power calculations. Power is currently only
calculated for time series of log-normally distributed data, which is
why the trend is expressed as a percentage.  
  

`recent_change`

A list with the following components:

- `n_year_fit` default = 5L

- `n_year_positive` default = 5L

A recent change will only be computed if there are at least `n_year_fit`
years of data in the recent period, of which at least `n_year_positive`
contain at least one non-censored measurement. This only affects
normally or log-normally distributed data.

The default values of 5L mirror the requirements for calculating the
change over the whole time series.

Note that in harsat versions \<= 1.0.2, `n_year_positive` was hard-wired
to 0L which occasionally led to pathological behaviour in the estimation
of the recent change. To replicate previous outputs as closely as
possible, set `n_year_positive` to 2L (the smallest value allowed to
avoid any pathologial behavour). The change is only likely to affect
long time series with infrequent sampling where most measurements in the
recent period are censored.
