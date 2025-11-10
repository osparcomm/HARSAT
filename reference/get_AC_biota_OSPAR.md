# Gets OSPAR threshold values for biota

Extends default extractor function get_AC\$biota for the few cases which
can not be handled in a straightforward manner. This mostly relates to
thresholds which are only applied when the typical species / lipid
contend is 'high' (currently defined as \>= 3%)

## Usage

``` r
get_AC_biota_OSPAR(data, AC, rt, export_all = FALSE)
```

## Arguments

- data:

  the data

- AC:

  a list of assessment criteria

- rt:

  the reference tables

- export_all:

  logical, default `FALSE`, if `TRUE` returns all data, otherwise, just
  the assessment criteria

## Value

by default, the assessment criteria, if `export_all` is set `TRUE`
returns all the data instead
