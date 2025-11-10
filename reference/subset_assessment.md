# Subsets an assessment object

Selects specific time series and simplifies the data, stations and
assessment components to match

## Usage

``` r
subset_assessment(assessment_obj, subset)
```

## Arguments

- assessment_obj:

  An assessment object resulting from a call to run_assessment.

- subset:

  A vector specifying the timeseries to be retained. An expression will
  be evaluated in the timeSeries component of assessment_obj; use
  'series' to identify individual timeseries.

## Value

a new assessment object, after applying the subset
