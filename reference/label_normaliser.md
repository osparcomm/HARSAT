# Text representation of normaliser information for plots and reports

Utilty function for providing inforamtion about normalistation to plots
in e.g. `plot_assessment` (where it is used with `label_unit`) or html
reports in e.g. `report_assessment`.

## Usage

``` r
label_normaliser(normaliser, value, unit, html = FALSE)
```

## Arguments

- normaliser:

  Character vector giving the names of the normaliser (ICES code rather
  than common name)

- value:

  Numeric vector giving the value of the normaliser

- unit:

  Character vector giving the unit of the normaliser

- html:

  Logical with TRUE returning an html representation for use in markdown
  and FALSE (default) returning an text expression for use in lattice
  (grid) graphics

## Value

A character string that can be used in markdown (`html = TRUE`) or used
in grid graphics ('html =
FALSE`). In the latter, need to turn the string into an expression by `parse(text
= result)\`.
