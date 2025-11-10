# Get the determinands to be assessed

Gets the determinands to be assessed from the determinand reference
table.

## Usage

``` r
ctsm_get_determinands(info)
```

## Arguments

- info:

  A list with at least the following two components:

  - `compartment` One of `"biota"`, `"sediment"` or `"water"`

  - `determinand` A data frame holding the determinand reference table

## Value

A character string containing the determinands to be assessed. The
function will fail with an error message if there are no such
determinands.

## Details

The determinands are taken from the column `biota_assess`,
`sediment_assess` or `water_assess` in the determinand reference table
(where the compartment is given by `info$compartment`). `TRUE` values
are determinands that are to be assessed.
