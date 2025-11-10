# Converts units, e.g., from mg/kg to ug/kg

Can accept non-standard units (e.g. for biological effects) provided
that `from` and `to` for these rows are identical (in which case no
attempt is made to convert.

## Usage

``` r
convert_units(conc, from, to)
```

## Arguments

- conc:

  the value to convert

- from:

  the current units

- to:

  the required units
