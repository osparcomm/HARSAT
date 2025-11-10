# Tidy the data

Reduces the size of the extraction by removing redundant variables. Any
ad-hoc changes will usually be made between
[`read_data`](http://osparcomm.github.io/HARSAT/reference/read_data.md)
and `tidy_data`. The output is in the correct format for
[`create_timeseries`](http://osparcomm.github.io/HARSAT/reference/create_timeseries.md).

## Usage

``` r
tidy_data(ctsm_obj)
```

## Arguments

- ctsm_obj:

  the time series data, typically returned from
  [`read_data`](http://osparcomm.github.io/HARSAT/reference/read_data.md).
