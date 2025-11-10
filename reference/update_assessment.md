# Update timeseries assessments

Refits models for particular timeseries, or fits new models when an
assessment is being done in chunks.

## Usage

``` r
update_assessment(ctsm_ob, subset = NULL, parallel = FALSE, ...)
```

## Arguments

- ctsm_ob:

  A HARSAT object resulting from a call to run_assessment

- subset:

  A vector specifying which timeseries assessements are to be updated or
  fit for the first time. An expression will be evaluated in the
  timeSeries component of ctsm_ob; use 'series' to identify individual
  timeseries.

- parallel:

  A logical which determines whether to use parallel computation;
  default = FALSE.

- ...:

  Extra arguments which are passed to assessment_engine. See details
  (which need to be written).
