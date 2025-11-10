# Reports the assessment of individual time series

Generates a series of html reports with, for each time series, meta
data, plots of the data with the fitted assessment model, statistical
summaries, and a simple interpretation of the fitted model.

## Usage

``` r
report_assessment(
  assessment_obj,
  subset = NULL,
  output_dir = ".",
  output_file = NULL,
  max_report = 100L
)
```

## Arguments

- assessment_obj:

  An assessment object resulting from a call to run_assessment

- subset:

  An optional vector specifying which timeseries are to be reported. An
  expression will be evaluated in the timeSeries component of
  assessment_obj; use 'series' to identify individual timeseries.

- output_dir:

  The output directory for the assessment plots (possibly supplied using
  'file.path'). The default is the working directory. The output
  directory must already exist.

- output_file:

  An alterntive file name to override the default. This is currently
  only implemented for a single report. If not supplied, the .html
  extension will be added.

- max_report:

  The maximum number of reports that will be generated. Defaults to 100.
  Each report is about 1MB in size and takes a few seconds to run, so
  this prevents a ridiculous number of reports being created.

## Value

A series of html files with, for each time series, meta data, plots of
the data with the fitted assessment model, statistical summaries, and a
simple interpretation of the fitted model.
