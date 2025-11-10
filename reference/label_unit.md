# Text representation of units for plots and reports

Utilty function for adding units to plots in e.g. `plot_assessment` or
html reports in e.g. `report_assessment`. Standard (recognised) units
are prettified, with non-standard units returned unchanged.

## Usage

``` r
label_unit(
  unit,
  basis = NA_character_,
  normaliser = NA_character_,
  normaliser_value = NA_real_,
  normaliser_unit = NA_character_,
  html = FALSE
)
```

## Arguments

- unit:

  Character vector giving the units

- basis:

  Character vector giving the bases. Must be one of "D", "L", "W" or
  NA_character\_

- normaliser:

  Character vector giving the names of the normaliser (if there are any)

- normaliser_value:

  Numeric vector giving the values of the normaliser (if there are any)

- normaliser_unit:

  Character vector giving the units of the normaliser (if there are any)

- html:

  Logical with TRUE returning an html representation for use in markdown
  and FALSE (default) returning an text expression for use in lattice
  (grid) graphics

## Value

A character vector that can be used in markdown (`html = TRUE`) or in
grid graphics ('html =
FALSE`). In the latter, need to turn the strings into an expression by `parse(text
= result)\`.
