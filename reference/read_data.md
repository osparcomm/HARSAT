# Read HARSAT data

Reads in contaminant and effects data, the station dictionary and
various reference tables. For data from the ICES webservice, it matches
data to stations in the station dictionary. It also allows the user to
set control parameters that dictate the assessment process.

## Usage

``` r
read_data(
  compartment = c("biota", "sediment", "water"),
  purpose = c("OSPAR", "HELCOM", "AMAP", "custom"),
  contaminants,
  stations,
  data_dir = ".",
  data_format = c("ICES", "external"),
  info_files = list(),
  info_dir = ".",
  extraction = NULL,
  max_year = NULL,
  oddity_dir = "oddities",
  control = list()
)
```

## Arguments

- compartment:

  A string: `"biota"`, `"sediment"` or `"water"`

- purpose:

  A string specifying whether to use the default set up for `"OSPAR"`,
  `"HELCOM"`, or `"AMAP"` or to use a customised setup `"custom"`

- contaminants:

  A file reference for the contaminant data

- stations:

  A file reference for the station data

- data_dir:

  The directory where the data files can be found (sometimes supplied
  using 'file.path'). Defaults to "."; i.e. the working directory.

- data_format:

  A string specifying whether the data were extracted from the ICES
  webservice (`"ICES"` - the default) or are in the simplified format
  designed for other data sources (`"external"`).

- info_files:

  A list of files specifying reference tables which override the
  defaults. See examples.

- info_dir:

  The directory where the reference tables can be found (sometimes
  supplied using 'file.path'). Defaults to "."; i.e. the working
  directory

- extraction:

  A date saying when the extraction was made. Optional. This should be
  provided according to ISO 8601; for example, 29 February 2024 should
  be supplied as "2024-02-29". If the contaminant data were extracted
  from the ICES webservice and the download file name has not been
  changed, the extraction data will be taken from the contaminant file
  name.

- max_year:

  An integer giving the last monitoring year that should be included in
  the assessment. Data from monitoring years after `max_year` will be
  deleted. If not specified `max_year` is taken to be the last
  monitoring year in the contaminant data file.

- oddity_dir:

  The directory where the 'oddities' will be written (sometimes supplied
  using 'file.path'). This directory (and subdirectories) will be
  created if it does not already exist.

- control:

  A list of control parameters that override the default values used to
  run the assessment. These include the reporting window; the way in
  which data are matched to stations following an ICES extraction;
  information about reporting regions, and so on. See Details.

## Value

A list with the following components:

- `call` The function call.

- `info` A list containing the reference tables and the control
  parameters.

- `data` A data frame containing the contaminant (and effects) data. For
  `external` data, this is identical to the input data file apart from
  some extra empty columns which have been added. For `ICES` data, some
  existing columns have been renamed (otherwise they are untouched) and
  some additional columns have been constructed. The key ones of these
  are:

  - `station_code` the code of the station in the station dictionary
    that best matches the data

  - `station_name` the name of the station

  - `species` (biota) the species based on `worms_accepted_name` where
    available and `speci_name` otherwise

  - `filtration` (water) whether the sample was `filtered` or
    `unfiltered` based on `method_pretreatment`

  - `retain` a logical indicating whether each record would have been
    retained under the previous ICES extraction protocol. For example,
    `retain` will be `FALSE` if the vflag entry is `"S"` or suspect.
    Records for which `retain == FALSE` are deleted later in `tidy_data`

- `stations`

## Details

### Control parameters

Many aspects of the assessment process can be controlled using
parameters which are stored in the `info` component of the harsat data
object. The default control values can be overwritten using the
`control` argument.

- `reporting_window` A scalar (default 6) which determines whether
  timeseries are excluded because they have no 'recent' data. Formally,
  timeseries are excluded if they have no data in the period
  `max_year - reporting_window + 1` and `max_year`, so the default
  approach is to exclude timeseries if they have no dat in the most
  recent six monitoring years. The value of 6 is chosen to match with
  Marine Strategy Framework Directive reporting periods.

- `region`

- `add_stations`

- `bivalve_spawning_season`

- `use_stage`

- `relative_uncertainty`

- `auxiliary` A list which allows flexibility in the treatment of
  auxiliary variables. At present, there is just one component
  `by_matrix`, a character vector that determines which auxiliary
  variables are matched to the contaminant data by `sample` and `matrix`
  as opposed to just `sample`. For sediment and water, the default is
  `all`; i.e. all variables are matched by `sample` and `matrix`. This
  ensures, for example, that sediment normalisers such as aluminium and
  organic carbon content are matched to chemical measurements in the
  same grain fraction. For biota, the default is
  `c("DRYWT%", "LIPIDWT%)`, so these variables are matched by `sample`
  and `matrix` and all other variables (e.g. LNMEA or %FEMALEPOP) are
  matched by `sample`. Thus, dry weight and lipid weight contents are
  matched to chemical measurements in the same tissue. However, mean
  length (which is usually the lenght of the whole organism) is matched
  to all tissue types.

### External data

If `data_format = "external"`, a simplified data and station file can be
supplied. See
[`vignette("external-file-format")`](http://osparcomm.github.io/HARSAT/articles/external-file-format.md)
for details.
