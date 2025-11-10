# Pretty representation of units for plots and reports

Utilty function for adding units to plots in e.g. `plot_assessment` or
html reports in e.g. `report_assessment`. Standard (recognised) units
are prettified, with non-standard units returned unchanged.

## Usage

``` r
pretty_unit(unit, html = FALSE)
```

## Arguments

- unit:

- html:

  A logical with TRUE returning an html representation for use in
  markdown and FALSE (default) returning an text expression for use in
  lattice (grid) graphics

## Value

A vector of character strings that can be used in markdown
(`html = TRUE`) or in grid graphics ('html =
FALSE`). In the latter, need to turn the strings into expressions by `parse(text
= result)\`.
