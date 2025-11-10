# Graphical summaries of an assessment

Generates a series of assessment plots for each time series. The plots
are exported as either png or pdf files.

## Usage

``` r
plot_assessment(
  assessment_obj,
  subset = NULL,
  output_dir = ".",
  file_type = c("data", "index", "auxiliary"),
  file_format = c("png", "pdf"),
  auxiliary = "default"
)
```

## Arguments

- assessment_obj:

  An assessment object resulting from a call to run_assessment

- subset:

  An optional vector specifying which timeseries are to be plotted. An
  expression will be evaluated in the timeSeries component of
  assessment_obj; use `series` to identify individual timeseries.

- output_dir:

  The output directory for the assessment plots (possibly supplied using
  `file.path`). The default is the working directory. The output
  directory must already exist.

- file_type:

  A character vector specifying the types of assessment plot. The
  default `c("data", "index", "auxiliary")` produces three plots for
  each time series. See details

- file_format:

  A character string specifying Whether the files should be png (the
  default) or pdf.

- auxiliary:

  A character string specifying the auxiliary variables plotted if
  `file_type = "auxiliary"`. See details

## Value

A series of png or pdf files with graphical summaries of an assessment.

## Details

### Types of assessment plots

- `file_type = "data"` shows the raw data with the fitted trend and
  pointwise two-sided 90% confidence bands

- `file_type = "index"` shows annual indices that summarise the data for
  each year with the fitted trend and pointwise two-sided 90% confidence
  bands

- `file_type = "auxiliary"` shows the raw data and key auxiliary
  variables see below)

### Auxiliary variables

The default (`auxiliary = "default"`) is to plot the following
variables:

- biota: determinand concentration, LNMEA (mean length), DRYWT% (dry
  weight content), LIPIDWT% (lipi weight content)

- sediment: non-normalised determinand concentration, normalised
  determinand concentration, AL (aluminium concentration), CORG (organic
  carbon content)

- water: no plots are generated at present

For biota, the determinand concentration will always be plotted, but it
is possible to change the three auxiliary variables. For example, to
plot WTMEA (mean weight) instead of LIPIDWT% you would set
`auxiliary = c("LNMEA", "WTMEA", "DRYWT%)`. For this to work, WTMEA must
previously have been specified as an auxiliary variable for the
determinand in question using the `biota_auxliary` column in the
determinand reference table. At present, there must always be three
auxiliary variables for biota.

For sediment, the non-normalised determinand concentration and the
normalised determinand concentration will always be plotted, but it is
possible to change the two auxiliary variables. For example, for metals
in sediment, you might set `auxiliary = c("AL", "LI")` to plot aluminium
and lithium concentrations instead of aluminium and organic carbon
concentrations. Again, for this to work, LI must previously have been
specified as an auxiliary variable for the determinand in question using
the `sediment_auxliary` column in the determinand reference table. At
present, there must always be two auxiliary variables for sediment.

At present, plots for only a limited range of auxiliary variables are
supported. More flexibility in these plots, such as changing the number
of auxiliary variables, is desirable and will emerge in due course.
