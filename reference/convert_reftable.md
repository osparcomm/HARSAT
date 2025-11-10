# Add OSPAR_subregion to simplified sediment threshold reference table

This is a utility function to expand a simplified OSPAR sediment
threshold table (which is much easier for the user to edit) into the
form required by `harsat`

## Usage

``` r
convert_reftable(input_file, output_file, export = TRUE)
```

## Arguments

- input_file:

  the input reference file

- output_file:

  the expanded reference file

- export:

  a boolean flag, if `FALSE`, the data is returned rather than being
  written to `output_file`

## Value

if `export` is `FALSE` (the default), returns the expanded data
