# Make key for annotating assessment plots

Provides the meta-data that accompanies each plot in `plot_assessment`:
compartment, station, units and data extraction. It also provides
extended information for plots involving multiple timeseries e.g. in
`plot_multidata`, `plot_multiassessment` and `plot_ratio`.

## Usage

``` r
make_key(series, info, type)
```

## Arguments

- series:

  List describing the timeseries in the plot. If there are multiple
  timeseries, then some entries of the list require multiple entries,
  one for each timeseries. See details.

- info:

  The standard `info` element of the harsat object.

- type:

  Character scalar (default `"data"`) describing the type of plot: other
  options are `"index"`, `"auxiliary"`, `"data_splom"`, `"index_mp"`,
  `"ratio_mp"`.

## Value

a list with four text (or expression) elements:

- `media`: compartment , species and matrix (biota), matrix (sediment),
  filtration (water), subseries

- `station`: station_name and station_longname

- `units`: an expression giving the units including any normalisation

- `extraction`: the extraction date

Both `media` and `units` can be complicated if there are multiple
timeseries, particularly if e.g. there are different normalisers
(sediment) or matrices (biological effects) involved \`

## Details

`series` is a list with giving information about a single timeseries
(when called by `plot_assessment`) or multiple related timeseries (when
called by `plot_multidata`, `plot_multiassessment`, and `plot_ratio`).
It must contain the following (compartment-specific) elements:

- `determinand` (all): a character vector giving the determinand for
  each timeseries; if there are multiple timeseries, then the vector
  must have names which are the corresponding series identifiers

- `matrix` (biota, sediment): a character vector giving the matrix for
  each timeseries; if there are multiple timeseries, then the vector
  must have names which are the corresponding series identifiers; the
  matrices might vary across timeseries, but this is currently likely
  only for biological effects

- `filtration` (water): a character vector giving the filtration method
  for each timeseries; if there are multiple timeseries, then the vector
  must have names which are the corresponding series identifiers; the
  filtration methods might vary across timeseries, but this is currently
  unlikely in practice

- `basis` (all): a character vector giving the basis for each
  timeseries; if there are multiple timeseries, then the vector must
  have names which are the corresponding series identifiers; the bases
  might vary across timeseries, but this is currently likely only for
  biological effects

- `species` (biota): a character giving the species (common to all
  timeseries)

- `subseries` (all): a character giving the subseries (common to all
  timeseries); can be `NA_character_`

In addition, if `info$normalise` is `TRUE`, then it must also contain:

- `normaliser`, `normaliser_value` and `normaliser_unit`: character,
  numeric, and character vectors giving the normaliser, it's value and
  the unit on which it is expressed for each timeseries; missing values
  are allowed if there is no normalisation for the timeseries involved
