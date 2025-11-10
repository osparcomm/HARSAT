# Extracts timeSeries

Gets the `timeSeries` component of a `harsat` object, optionally adding
extra information about each station

## Usage

``` r
get_timeseries(harsat_obj, add = TRUE)
```

## Arguments

- harsat_obj:

  A `harsat` object following a call to
  [`create_timeseries`](http://osparcomm.github.io/HARSAT/reference/create_timeseries.md)

- add:

  logical (default `TRUE`); `TRUE` adds the `station_name` and `country`
  associated with each station; `FALSE` simply returns the existing
  timeseries information

## Value

A data.frame containing the `timeSeries` component with (optionally)
extra information about each station. The `series` column is the
identifier of each timeseries.
