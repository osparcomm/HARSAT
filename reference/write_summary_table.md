# Write assessment summary to a csv file

Creates a data frame summarising the assessment of each time series and
writes it to a csv file. The summary includes:

- meta-data such as the monitoring location and number of years of data
  for each time series

- the fitted values in the last monitoring year with associated upper
  one-sided 95% confidence limits

- the trend assessments (p-values and trend estimates)

- the status assessments (if there any thresholds)

- (optionally) a symbology summarising the trend (shape) and status
  (colour) of each time series. This is experimental.

## Usage

``` r
write_summary_table(
  harsat_obj,
  output_file = NULL,
  output_dir = ".",
  export = TRUE,
  threshold_groups = NULL,
  collapse_AC = lifecycle::deprecated(),
  extra_output = NULL,
  symbology = NULL,
  symbology_control = list(),
  determinandGroups = NULL,
  append = FALSE
)
```

## Arguments

- harsat_obj:

  A harsat object following a call to `run_assessment`.

- output_file:

  The name of the output csv file. If using NULL, the file will be
  called `biota_summary.csv`, `sediment_summary.csv` or
  `water_summary.csv` as appropriate. By default the file will be
  written to the working directory. If a file name is provided, a path
  to the output file can also be provided (e.g. using `file.path`). The
  \`output_dir“ option can also be used to specify the output file
  directory.

- output_dir:

  The output directory for `output_file`. The default is the working
  directory. Any file path provided in `output_file`, will be appended
  to `output_dir`. The resulting output directory must already exist.

- export:

  Logical. `TRUE` (the default) writes the summary table to a csv file.
  `FALSE` returns the summary table as an R object (and does not write
  to a csv file).

- threshold_groups:

  A names list of valid thresholds that allows thresholds of the same
  'type' to be reported together. See details.

- collapse_AC:

  **\[deprecated\]** Use `threshold_groups` instead.

- extra_output:

  A character vector specifying extra summary metrics to be included in
  the output. Currently only recognises "power" to give the seven power
  metrics computed for lognormally distributed data. Defaults to `NULL`;
  i.e. no extra output.

- symbology:

  Experimental. A character string "default" or a user-defined function
  that specifies a symbology typically used to characterise the patterns
  of change in and the status of each time series. Defaults to `NULL`;
  i.e. no symbology. Multiple symbologies can be applied. See details.

- symbology_control:

  Experimental. A named list of control options for the symbology. See
  details.

- determinandGroups:

  optional, a list specifying `labels` and `levels` to rename the
  existing determinand groups. The life of this argument is limited.

- append:

  Logical. `FALSE` (the default) overwrites any existing summary file.
  `TRUE` appends data to it, creating it if it does not yet exist.

## Value

a summary object, when `export` is `FALSE`

## Default symbology

`symbology = "default"` calls a pre-defined symbology that generates a
'shape' and a 'colour' to characterise the status of each time series.
Its behaviour is controlled using `symbology_control`, a named list with
the following elements:

- `shape`: a list with names `none`, `mean`, `flat`, `up`, `down` giving
  the shape associated with each pattern of change. Here, `none`
  corresponds to insufficient data to fit a parametric model; `mean` to
  sufficient data to fit a parametric model but not to assess for
  trends; `flat` to no significant change in level (concentration) over
  time; `up` to a \#' significant increase in level over time; `down` to
  a significant decrease in level over time. Their default values are
  `"small_open_circle"`, `"small_filled_circle"`,
  `"large_filled_circle"`, `"upward_triangle"`, `"downward_triangle"`

- `alpha` is the size of the test for change; default = `0.05`

- `change` determines whether the change is based on the recent time
  window (typically the last twenty years) or the whole time series;
  options `"recent"` (default) and `"overall"`

- `colour` (default `NULL`) is a named list that characterises the
  status of a time series based on specified thresholds. The list names
  must match (a subset of) the names of the thresholds used in the
  assessment or, if the thresholds have been grouped, the group names.
  Each threshold must have two elements: `below` gives the colour if the
  time series is significantly below the threshold (p \< 0.05); `above`
  gives the colour otherwise. If multiple thresholds are used, they must
  be ordered from best to worst status. See examples

- `no_threshold` (default `"black"`) is the colour used when no
  thresholds are applied to a time series. Another option might be to
  use `NA_character_`

- `adjust_nonparam` (default `TRUE`) is a logical that allows the
  symbology to be adjusted for short time series (often dominated by
  less-than values) where a non-parametric test for status can be
  applied

- `names` (default `list(colour = "colour", shape = "shape")`) allows
  the names of the symbology columns in the summary table to be
  adjusted; this can be important if multiple symbologies are applied

## Custom symbologies

Users can apply custom symbologies by letting `symbology` be a
user-supplied function of the form `fn(summary, info, control)` where:

- `summary` is the summary table before applying the symbology; for
  convenience, there is an additional column `method` which doesn't
  appear in the final summary table but takes, in particular, values
  `"none"` and `"mean"` corresponding respectively to no parametric
  model and insufficient data to fit a trend (see `shape` above)

- `info` contains the contents of `harsat_obj$info`; i.e. all the
  reference tables and additional information about the assessment. For
  convenience, it also contains a temporary element `.threshold_group`
  which has the names of the threshold groups (which can differ from the
  thresholds themselves)

- `control` contains `symbology_control` and allows the user to pass
  additional information to the function

The output of the function must be a data frame with one column called
`series` that contains the series identifier of each time series (not
necessarily in the same order as the summary table) and symbology
columns, but they must not share any of the existing names in the
summary table.

See the examples for more inspiration.

## Multiple symbologies

Multiple symbologies can be applied by specifying a named list whch can
be a mixture of default symbologies and custom symbologies.
`symbology_control` must then be a named list (with the same names)
giving control information for each symbology. It is important to ensure
that each symbology gives output columns with different names. See
examples.

## Examples

``` r

# Default symbology with one threshold: the EQS. The colour will be "green"
# if the time series is significantly below the EQS in the last monitoring
# year and "red" otherwise 
if (FALSE) { # \dontrun{
write_summary_table(
  water_assessment,
  symbology = "default",
  symbology_control = list(
    colour = list(EQS = list(below = "green", above = "red"))
  )
)
} # }

# Now applied using the overall change instead of the recent change.
if (FALSE) { # \dontrun{
write_summary_table(
  water_assessment,
  symbology = "default",
  symbology_control = list(
    colour = list(EQS = list(below = "green", above = "red")),
    change = "overall"
  )
)
} # }

# If we only want to change one shape, then we only need to specify that one
if (FALSE) { # \dontrun{
write_summary_table(
  water_assessment,
  symbology = "default",
  symbology_control = list(
    colour = list(EQS = list(below = "green", above = "red")),
    shape = list(flat = "square")
  )
)
} # }

# Assessment thresholds grouped into BAC and EAC equivalents.
# Symbology now has two thresholds giving:
# "blue" if significantly below the BAC
# "orange" if not significantly below the BAC and there is no EAC
# "green" if significantly below the EAC but not the BAC
# "red" otherwise  
if (FALSE) { # \dontrun{
write_summary_table(
  sediment_assessment,
  threshold_groups = list(
    BAC = "BAC", 
    EAC = c("EAC", "ERL", "EQS", "FEQG")
  ),
  symbology = "default", 
  symbology_control = list(
    colour = list(
      BAC = list(below = "blue", above = "orange"),
      EAC = list(below = "green", above = "red")
    ), 
  )
)
} # }

# Assessment thresholds grouped into BAC and EAC equivalents. Human health
# thresholds grouped as HQS
# Two symbologies applied, one for environmental thresholds, the other for
# health thresholds
# Note the named lists and that the output names are specified 
if (FALSE) { # \dontrun{
write_summary_table(
  biota_assessment,
  threshold_groups = list(
    BAC = c("BAC", "NRC"),
    EAC = c("EAC", "FEQG", "LRC", "QSsp"), 
    HQS = c("MPC", "QShh")
  ),
  symbology = list(env = "default", health = "default"), 
  symbology_control = list(
    env = list(
      colour = list(
        BAC = list(below = "blue", above = "orange"),
        EAC = list(below = "green", above = "red")
      ), 
      names = list(shape = "shape_env", colour = "colour_env")
    ),
    health = list(
      colour = list(HQS = list(below = "green", above = "red")), 
      names = list(shape = "shape_health", colour = "colour_health")
    )
  )
)
} # }

# Custom symbology that only reports time series where there is sufficient
# information to assess trends and which colours the time series by whether, 
# for each determinand, mean concentrations in the last monitoring year are 
# below or above the median concentration observed across time series 
if (FALSE) { # \dontrun{
symbology_user <- function(summary, info, control) {
  summary <- dplyr::mutate(
    summary,
    shape = dplyr::case_when(
      is.na(p_overall_change)   ~ NA_character_,
      p_overall_change > 0.05   ~ "circle",
      overall_change > 0        ~ "upward_triangle",
      overall_change < 0        ~ "downward_triangle"
    )
  )
  summary <- summary |> 
    dplyr::group_by(determinand) |>
    dplyr::mutate(
      .shape = !is.na(shape),
      colour = dplyr::case_when(
        !.shape                                          ~ NA_character_,
        mean_last_year <= median(mean_last_year[.shape]) ~ "blue",
       mean_last_year > median(mean_last_year[.shape])  ~ "red"
      )
    ) |>
    dplyr::ungroup()
  summary[c("series", "shape", "colour")]
}

write_summary_table(
  biota_assessment,
  symbology = symbology_user
)
} # }

```
